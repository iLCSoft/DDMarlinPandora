/**
 *  @file   MarlinPandora/src/TrackCreator.cc
 * 
 *  @brief  Implementation of the track creator class.
 * 
 *  $Log: $
 */

#include "marlin/Global.h"
#include "marlin/Processor.h"

#include "EVENT/LCCollection.h"
#include "EVENT/ReconstructedParticle.h"
#include "EVENT/Vertex.h"
#include "UTIL/ILDConf.h"

#include "gear/BField.h"
#include "gear/CalorimeterParameters.h"
#include "gear/PadRowLayout2D.h"
#include "gear/TPCParameters.h"
#include "gear/FTDParameters.h"
#include "gear/FTDLayerLayout.h"

#include "DDPandoraPFANewProcessor.h"
#include "TrackCreator.h"
#include "Pandora/PdgTable.h"

#include <algorithm>
#include <cmath>
#include <limits>

TrackCreator::TrackCreator(const Settings &settings, const pandora::Pandora *const pPandora) :
    m_settings(settings),
    m_pPandora(pPandora),
    m_bField(marlin::Global::GEAR->getBField().at(gear::Vector3D(0., 0., 0.)).z()),
    m_tpcInnerR(marlin::Global::GEAR->getTPCParameters().getPadLayout().getPlaneExtent()[0]),
    m_tpcOuterR(marlin::Global::GEAR->getTPCParameters().getPadLayout().getPlaneExtent()[1]),
    m_tpcMaxRow(marlin::Global::GEAR->getTPCParameters().getPadLayout().getNRows()),
    m_tpcZmax(marlin::Global::GEAR->getTPCParameters().getMaxDriftLength()),
    m_eCalBarrelInnerSymmetry(marlin::Global::GEAR->getEcalBarrelParameters().getSymmetryOrder()),
    m_eCalBarrelInnerPhi0(marlin::Global::GEAR->getEcalBarrelParameters().getPhi0()),
    m_eCalBarrelInnerR(marlin::Global::GEAR->getEcalBarrelParameters().getExtent()[0]),
    m_eCalEndCapInnerZ(marlin::Global::GEAR->getEcalEndcapParameters().getExtent()[2])
{
    // fg: FTD description in GEAR has changed ...
    try
    {
        m_ftdInnerRadii = marlin::Global::GEAR->getGearParameters("FTD").getDoubleVals("FTDInnerRadius");
        m_ftdOuterRadii = marlin::Global::GEAR->getGearParameters("FTD").getDoubleVals("FTDOuterRadius");
        m_ftdZPositions = marlin::Global::GEAR->getGearParameters("FTD").getDoubleVals("FTDZCoordinate");
        m_nFtdLayers = m_ftdZPositions.size();
    }
    catch (gear::UnknownParameterException &)
    {
        const gear::FTDLayerLayout &ftdLayerLayout(marlin::Global::GEAR->getFTDParameters().getFTDLayerLayout());
        streamlog_out( DEBUG2 ) << " Filling FTD parameters from gear::FTDParameters - n layers: " << ftdLayerLayout.getNLayers() << std::endl;

        for(unsigned int i = 0, N = ftdLayerLayout.getNLayers(); i < N; ++i)
        {
            // Create a disk to represent even number petals front side
            m_ftdInnerRadii.push_back(ftdLayerLayout.getSensitiveRinner(i));
            m_ftdOuterRadii.push_back(ftdLayerLayout.getMaxRadius(i));

            // Take the mean z position of the staggered petals
            const double zpos(ftdLayerLayout.getZposition(i));
            m_ftdZPositions.push_back(zpos);
            streamlog_out( DEBUG2 ) << "     layer " << i << " - mean z position = " << zpos << std::endl;
        }

        m_nFtdLayers = m_ftdZPositions.size() ;
    }

    // Check tpc parameters
    if ((std::fabs(m_tpcZmax) < std::numeric_limits<float>::epsilon()) || (std::fabs(m_tpcInnerR) < std::numeric_limits<float>::epsilon())
        || (std::fabs(m_tpcOuterR - m_tpcInnerR) < std::numeric_limits<float>::epsilon()))
    {
        throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);
    }

    m_cosTpc = m_tpcZmax / std::sqrt(m_tpcZmax * m_tpcZmax + m_tpcInnerR * m_tpcInnerR);

    // Check ftd parameters
    if ((0 == m_nFtdLayers) || (m_nFtdLayers != m_ftdInnerRadii.size()) || (m_nFtdLayers != m_ftdOuterRadii.size()))
    {
        throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);
    }

    for (unsigned int iFtdLayer = 0; iFtdLayer < m_nFtdLayers; ++iFtdLayer)
    {
        if ((std::fabs(m_ftdOuterRadii[iFtdLayer]) < std::numeric_limits<float>::epsilon()) ||
            (std::fabs(m_ftdInnerRadii[iFtdLayer]) < std::numeric_limits<float>::epsilon()))
        {
            throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);
        }
    }

    m_tanLambdaFtd = m_ftdZPositions[0] / m_ftdOuterRadii[0];

    // Calculate etd and set parameters
    // fg: make SET and ETD optional - as they might not be in the model ...
    try
    {
        const DoubleVector &etdZPositions(marlin::Global::GEAR->getGearParameters("ETD").getDoubleVals("ETDLayerZ"));
        const DoubleVector &setInnerRadii(marlin::Global::GEAR->getGearParameters("SET").getDoubleVals("SETLayerRadius"));

        if (etdZPositions.empty() || setInnerRadii.empty())
            throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);

        m_minEtdZPosition = *(std::min_element(etdZPositions.begin(), etdZPositions.end()));
        m_minSetRadius = *(std::min_element(setInnerRadii.begin(), setInnerRadii.end()));
    }
    catch(gear::UnknownParameterException &)
    {
        streamlog_out(WARNING) << " ETDLayerZ or SETLayerRadius parameters missing from GEAR parameters!" << std::endl
                               << "     -> both will be set to " << std::numeric_limits<float>::quiet_NaN() << std::endl;

        //fg: Set them to NAN, so that they cannot be used to set   trackParameters.m_reachesCalorimeter = true;
        m_minEtdZPosition = std::numeric_limits<float>::quiet_NaN();
        m_minSetRadius = std::numeric_limits<float>::quiet_NaN();
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

TrackCreator::~TrackCreator()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode TrackCreator::CreateTrackAssociations(const EVENT::LCEvent *const pLCEvent)
{
    PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, this->ExtractKinks(pLCEvent));
    PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, this->ExtractProngsAndSplits(pLCEvent));
    PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, this->ExtractV0s(pLCEvent));

    return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode TrackCreator::ExtractKinks(const EVENT::LCEvent *const pLCEvent)
{
    for (StringVector::const_iterator iter = m_settings.m_kinkVertexCollections.begin(), iterEnd = m_settings.m_kinkVertexCollections.end();
        iter != iterEnd; ++iter)
    {
        try
        {
            const EVENT::LCCollection *pKinkCollection = pLCEvent->getCollection(*iter);

            for (int i = 0, iMax = pKinkCollection->getNumberOfElements(); i < iMax; ++i)
            {
                try
                {
                    EVENT::Vertex *pVertex = dynamic_cast<Vertex*>(pKinkCollection->getElementAt(i));

                    if (NULL == pVertex)
                        throw EVENT::Exception("Collection type mismatch");

                    EVENT::ReconstructedParticle *pReconstructedParticle = pVertex->getAssociatedParticle();
                    const EVENT::TrackVec &trackVec(pReconstructedParticle->getTracks());

                    if (this->IsConflictingRelationship(trackVec))
                        continue;

                    const int vertexPdgCode(pReconstructedParticle->getType());

                    // Extract the kink vertex information
                    for (unsigned int iTrack = 0, nTracks = trackVec.size(); iTrack < nTracks; ++iTrack)
                    {
                        EVENT::Track *pTrack = trackVec[iTrack];
                        (0 == iTrack) ? m_parentTrackList.insert(pTrack) : m_daughterTrackList.insert(pTrack);
                        streamlog_out(DEBUG) << "KinkTrack " << iTrack << ", nHits " << pTrack->getTrackerHits().size() << std::endl;

                        int trackPdgCode = pandora::UNKNOWN_PARTICLE_TYPE;

                        if (0 == iTrack)
                        {
                            trackPdgCode = vertexPdgCode;
                        }
                        else
                        {
                            switch (vertexPdgCode)
                            {
                            case pandora::PI_PLUS :
                            case pandora::K_PLUS :
                                trackPdgCode = pandora::MU_PLUS;
                                break;
                            case pandora::PI_MINUS :
                            case pandora::K_MINUS :
                                trackPdgCode = pandora::MU_MINUS;
                                break;
                            case pandora::HYPERON_MINUS_BAR :
                            case pandora::SIGMA_PLUS :
                                trackPdgCode = pandora::PI_PLUS;
                                break;
                            case pandora::SIGMA_MINUS :
                            case pandora::HYPERON_MINUS :
                                trackPdgCode = pandora::PI_PLUS;
                                break;
                            default :
                                (pTrack->getOmega() > 0) ? trackPdgCode = pandora::PI_PLUS : trackPdgCode = pandora::PI_MINUS;
                                break;
                            }
                        }

                        m_trackToPidMap.insert(TrackToPidMap::value_type(pTrack, trackPdgCode));

                        if (0 == m_settings.m_shouldFormTrackRelationships)
                            continue;

                        // Make track parent-daughter relationships
                        if (0 == iTrack)
                        {
                            for (unsigned int jTrack = iTrack + 1; jTrack < nTracks; ++jTrack)
                            {
                                PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::SetTrackParentDaughterRelationship(*m_pPandora,
                                    pTrack, trackVec[jTrack]));
                            }
                        }

                        // Make track sibling relationships
                        else
                        {
                            for (unsigned int jTrack = iTrack + 1; jTrack < nTracks; ++jTrack)
                            {
                                PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::SetTrackSiblingRelationship(*m_pPandora,
                                    pTrack, trackVec[jTrack]));
                            }
                        }
                    }
                }
                catch (EVENT::Exception &exception)
                {
                    streamlog_out(WARNING) << "Failed to extract kink vertex: " << exception.what() << std::endl;
                }
            }
        }
        catch (EVENT::Exception &exception)
        {
            streamlog_out(DEBUG5) << "Failed to extract kink vertex collection: " << *iter << ", " << exception.what() << std::endl;
        }
    }

    return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode TrackCreator::ExtractProngsAndSplits(const EVENT::LCEvent *const pLCEvent)
{
    for (StringVector::const_iterator iter = m_settings.m_prongSplitVertexCollections.begin(), iterEnd = m_settings.m_prongSplitVertexCollections.end();
        iter != iterEnd; ++iter)
    {
        try
        {
            const EVENT::LCCollection *pProngOrSplitCollection = pLCEvent->getCollection(*iter);

            for (int i = 0, iMax = pProngOrSplitCollection->getNumberOfElements(); i < iMax; ++i)
            {
                try
                {
                    EVENT::Vertex *pVertex = dynamic_cast<Vertex*>(pProngOrSplitCollection->getElementAt(i));

                    if (NULL == pVertex)
                        throw EVENT::Exception("Collection type mismatch");

                    EVENT::ReconstructedParticle *pReconstructedParticle = pVertex->getAssociatedParticle();
                    const EVENT::TrackVec &trackVec(pReconstructedParticle->getTracks());

                    if (this->IsConflictingRelationship(trackVec))
                        continue;

                    // Extract the prong/split vertex information
                    for (unsigned int iTrack = 0, nTracks = trackVec.size(); iTrack < nTracks; ++iTrack)
                    {
                        EVENT::Track *pTrack = trackVec[iTrack];
                        (0 == iTrack) ? m_parentTrackList.insert(pTrack) : m_daughterTrackList.insert(pTrack);
                        streamlog_out(DEBUG) << "Prong or Split Track " << iTrack << ", nHits " << pTrack->getTrackerHits().size() << std::endl;

                        if (0 == m_settings.m_shouldFormTrackRelationships)
                            continue;

                        // Make track parent-daughter relationships
                        if (0 == iTrack)
                        {
                            for (unsigned int jTrack = iTrack + 1; jTrack < nTracks; ++jTrack)
                            {
                                PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::SetTrackParentDaughterRelationship(*m_pPandora,
                                    pTrack, trackVec[jTrack]));
                            }
                        }

                        // Make track sibling relationships
                        else
                        {
                            for (unsigned int jTrack = iTrack + 1; jTrack < nTracks; ++jTrack)
                            {
                                PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::SetTrackSiblingRelationship(*m_pPandora,
                                    pTrack, trackVec[jTrack]));
                            }
                        }
                    }
                }
                catch (EVENT::Exception &exception)
                {
                    streamlog_out(WARNING) << "Failed to extract prong/split vertex: " << exception.what() << std::endl;
                }
            }
        }
        catch (EVENT::Exception &exception)
        {
            streamlog_out(DEBUG5) << "Failed to extract prong/split vertex collection: " << *iter << ", " << exception.what() << std::endl;
        }
    }

    return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode TrackCreator::ExtractV0s(const EVENT::LCEvent *const pLCEvent)
{
    for (StringVector::const_iterator iter = m_settings.m_v0VertexCollections.begin(), iterEnd = m_settings.m_v0VertexCollections.end();
        iter != iterEnd; ++iter)
    {
        try
        {
            const EVENT::LCCollection *pV0Collection = pLCEvent->getCollection(*iter);

            for (int i = 0, iMax = pV0Collection->getNumberOfElements(); i < iMax; ++i)
            {
                try
                {
                    EVENT::Vertex *pVertex = dynamic_cast<Vertex*>(pV0Collection->getElementAt(i));

                    if (NULL == pVertex)
                        throw EVENT::Exception("Collection type mismatch");

                    EVENT::ReconstructedParticle *pReconstructedParticle = pVertex->getAssociatedParticle();
                    const EVENT::TrackVec &trackVec(pReconstructedParticle->getTracks());

                    if (this->IsConflictingRelationship(trackVec))
                        continue;

                    // Extract the v0 vertex information
                    const int vertexPdgCode(pReconstructedParticle->getType());

                    for (unsigned int iTrack = 0, nTracks = trackVec.size(); iTrack < nTracks; ++iTrack)
                    {
                        EVENT::Track *pTrack = trackVec[iTrack];
                        m_v0TrackList.insert(pTrack);
                        streamlog_out(DEBUG) << "V0Track " << iTrack << ", nHits " << pTrack->getTrackerHits().size() << std::endl;

                        int trackPdgCode = pandora::UNKNOWN_PARTICLE_TYPE;

                        switch (vertexPdgCode)
                        {
                        case pandora::PHOTON :
                            (pTrack->getOmega() > 0) ? trackPdgCode = pandora::E_PLUS : trackPdgCode = pandora::E_MINUS;
                            break;
                        case pandora::LAMBDA :
                            (pTrack->getOmega() > 0) ? trackPdgCode = pandora::PROTON : trackPdgCode = pandora::PI_MINUS;
                            break;
                        case pandora::LAMBDA_BAR :
                            (pTrack->getOmega() > 0) ? trackPdgCode = pandora::PI_PLUS : trackPdgCode = pandora::PROTON_BAR;
                            break;
                        case pandora::K_SHORT :
                            (pTrack->getOmega() > 0) ? trackPdgCode = pandora::PI_PLUS : trackPdgCode = pandora::PI_MINUS;
                            break;
                        default :
                            (pTrack->getOmega() > 0) ? trackPdgCode = pandora::PI_PLUS : trackPdgCode = pandora::PI_MINUS;
                            break;
                        }

                        m_trackToPidMap.insert(TrackToPidMap::value_type(pTrack, trackPdgCode));

                        if (0 == m_settings.m_shouldFormTrackRelationships)
                            continue;

                        // Make track sibling relationships
                        for (unsigned int jTrack = iTrack + 1; jTrack < nTracks; ++jTrack)
                        {
                            PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::SetTrackSiblingRelationship(*m_pPandora,
                                pTrack, trackVec[jTrack]));
                        }
                    }
                }
                catch (EVENT::Exception &exception)
                {
                    streamlog_out(WARNING) << "Failed to extract v0 vertex: " << exception.what() << std::endl;
                }
            }
        }
        catch (EVENT::Exception &exception)
        {
            streamlog_out(DEBUG5) << "Failed to extract v0 vertex collection: " << *iter << ", " << exception.what() << std::endl;
        }
    }

    return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TrackCreator::IsConflictingRelationship(const EVENT::TrackVec &trackVec) const
{
    for (unsigned int iTrack = 0, nTracks = trackVec.size(); iTrack < nTracks; ++iTrack)
    {
        EVENT::Track *pTrack = trackVec[iTrack];

        if (this->IsDaughter(pTrack) || this->IsParent(pTrack) || this->IsV0(pTrack))
            return true;
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode TrackCreator::CreateTracks(EVENT::LCEvent *pLCEvent)
{
    for (StringVector::const_iterator iter = m_settings.m_trackCollections.begin(), iterEnd = m_settings.m_trackCollections.end();
        iter != iterEnd; ++iter)
    {
        try
        {
            const EVENT::LCCollection *pTrackCollection = pLCEvent->getCollection(*iter);

            for (int i = 0, iMax = pTrackCollection->getNumberOfElements(); i < iMax; ++i)
            {
                try
                {
                    EVENT::Track *pTrack = dynamic_cast<Track*>(pTrackCollection->getElementAt(i));

                    if (NULL == pTrack)
                        throw EVENT::Exception("Collection type mismatch");

                    int minTrackHits = m_settings.m_minTrackHits;
                    const float tanLambda(std::fabs(pTrack->getTanLambda()));

                    if (tanLambda > m_tanLambdaFtd)
                    {
                        int expectedFtdHits(0);

                        for (unsigned int iFtdLayer = 0; iFtdLayer < m_nFtdLayers; ++iFtdLayer)
                        {
                            if ((tanLambda > m_ftdZPositions[iFtdLayer] / m_ftdOuterRadii[iFtdLayer]) &&
                                (tanLambda < m_ftdZPositions[iFtdLayer] / m_ftdInnerRadii[iFtdLayer]))
                            {
                                expectedFtdHits++;
                            }
                        }

                        minTrackHits = std::max(m_settings.m_minFtdTrackHits, expectedFtdHits);
                    }

                    const int nTrackHits(static_cast<int>(pTrack->getTrackerHits().size()));

                    if ((nTrackHits < minTrackHits) || (nTrackHits > m_settings.m_maxTrackHits))
                        continue;

                    // Proceed to create the pandora track
                    PandoraApi::Track::Parameters trackParameters;
                    trackParameters.m_d0 = pTrack->getD0();
                    trackParameters.m_z0 = pTrack->getZ0();
                    trackParameters.m_pParentAddress = pTrack;

                    // By default, assume tracks are charged pions
                    const float signedCurvature(pTrack->getOmega());
                    trackParameters.m_particleId = (signedCurvature > 0) ? pandora::PI_PLUS : pandora::PI_MINUS;
                    trackParameters.m_mass = pandora::PdgTable::GetParticleMass(pandora::PI_PLUS);

                    // Use particle id information from V0 and Kink finders
                    TrackToPidMap::const_iterator iter = m_trackToPidMap.find(pTrack);

                    if(iter != m_trackToPidMap.end())
                    {
                        trackParameters.m_particleId = (*iter).second;
                        trackParameters.m_mass = pandora::PdgTable::GetParticleMass((*iter).second);
                    }

                    if (0.f != signedCurvature)
                        trackParameters.m_charge = static_cast<int>(signedCurvature / std::fabs(signedCurvature));

                    this->GetTrackStates(pTrack, trackParameters);
                    this->TrackReachesECAL(pTrack, trackParameters);
                    this->DefineTrackPfoUsage(pTrack, trackParameters);

                    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::Track::Create(*m_pPandora, trackParameters));
                    m_trackVector.push_back(pTrack);
                }
                catch (pandora::StatusCodeException &statusCodeException)
                {
                    streamlog_out(ERROR) << "Failed to extract a track: " << statusCodeException.ToString() << std::endl;
                }
                catch (EVENT::Exception &exception)
                {
                    streamlog_out(WARNING) << "Failed to extract a vertex: " << exception.what() << std::endl;
                }
            }
        }
        catch (EVENT::Exception &exception)
        {
            streamlog_out(WARNING) << "Failed to extract track collection: " << *iter << ", " << exception.what() << std::endl;
        }
    }

    return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackCreator::GetTrackStates(const EVENT::Track *const pTrack, PandoraApi::Track::Parameters &trackParameters) const
{
    if (m_settings.m_useOldTrackStateCalculation > 0)
        return this->GetTrackStatesOld(pTrack, trackParameters);

    const TrackState *pTrackState = pTrack->getTrackState(TrackState::AtIP);

    if (!pTrackState)
        throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_INITIALIZED);

    const double pt(m_bField * 2.99792e-4 / std::fabs(pTrackState->getOmega()));
    trackParameters.m_momentumAtDca = pandora::CartesianVector(std::cos(pTrackState->getPhi()), std::sin(pTrackState->getPhi()), pTrackState->getTanLambda()) * pt;

    this->CopyTrackState(pTrack->getTrackState(TrackState::AtFirstHit), trackParameters.m_trackStateAtStart);

    //fg: curling TPC tracks have pointers to track segments stored -> need to get track states from last segment!
    const EVENT::Track *pEndTrack = (pTrack->getTracks().empty() ?  pTrack  :  pTrack->getTracks().back());

    this->CopyTrackState(pEndTrack->getTrackState(TrackState::AtLastHit), trackParameters.m_trackStateAtEnd);
    this->CopyTrackState(pEndTrack->getTrackState(TrackState::AtCalorimeter), trackParameters.m_trackStateAtCalorimeter);

    trackParameters.m_isProjectedToEndCap = ((std::fabs(trackParameters.m_trackStateAtCalorimeter.Get().GetPosition().GetZ()) < m_eCalEndCapInnerZ) ? false : true);

    // TODO minGenericTime * particleEnergy / 300.f;
    trackParameters.m_timeAtCalorimeter = 0.f;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackCreator::CopyTrackState(const TrackState *const pTrackState, pandora::InputTrackState &inputTrackState) const
{
    if (!pTrackState)
        throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_INITIALIZED);

    const double pt(m_bField * 2.99792e-4 / std::fabs(pTrackState->getOmega()));

    const double px(pt * std::cos(pTrackState->getPhi()));
    const double py(pt * std::sin(pTrackState->getPhi()));
    const double pz(pt * pTrackState->getTanLambda());

    const double xs(pTrackState->getReferencePoint()[0]);
    const double ys(pTrackState->getReferencePoint()[1]);
    const double zs(pTrackState->getReferencePoint()[2]);

    inputTrackState = pandora::TrackState(xs, ys, zs, px, py, pz);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackCreator::GetTrackStatesOld(const EVENT::Track *const pTrack, PandoraApi::Track::Parameters &trackParameters) const
{
    pandora::Helix *pHelixFit = new pandora::Helix(pTrack->getPhi(), pTrack->getD0(), pTrack->getZ0(), pTrack->getOmega(), pTrack->getTanLambda(), m_bField);
    trackParameters.m_momentumAtDca = pHelixFit->GetMomentum();

    const EVENT::TrackerHitVec &trackerHitvec(pTrack->getTrackerHits());
    float zMin(std::numeric_limits<float>::max()), zMax(-std::numeric_limits<float>::max());

    for (int iz = 0, nTrackHits = trackerHitvec.size(); iz < nTrackHits - 1; ++iz)
    {
        const float hitZ(trackerHitvec[iz]->getPosition()[2]);

        if (hitZ > zMax)
            zMax = hitZ;

        if (hitZ < zMin)
            zMin = hitZ;
    }

    const int signPz((pHelixFit->GetMomentum().GetZ() > 0.f) ? 1 : -1);
    const float zStart((signPz > 0) ? zMin : zMax);
    const float zEnd((signPz > 0) ? zMax : zMin);

    pandora::CartesianVector startPosition(0.f, 0.f, 0.f), startMomentum(0.f, 0.f, 0.f);
    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, pHelixFit->GetPointInZ(zStart, pHelixFit->GetReferencePoint(), startPosition));
    startMomentum = pHelixFit->GetExtrapolatedMomentum(startPosition);
    trackParameters.m_trackStateAtStart = pandora::TrackState(startPosition, startMomentum);

    pandora::CartesianVector endPosition(0.f, 0.f, 0.f), endMomentum(0.f, 0.f, 0.f);
    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, pHelixFit->GetPointInZ(zEnd, pHelixFit->GetReferencePoint(), endPosition));
    endMomentum = pHelixFit->GetExtrapolatedMomentum(endPosition);
    trackParameters.m_trackStateAtEnd = pandora::TrackState(endPosition, endMomentum);

    this->GetECalProjectionOld(pHelixFit, signPz, trackParameters);

    delete pHelixFit;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackCreator::GetECalProjectionOld(const pandora::Helix *const pHelix, const int signPz, PandoraApi::Track::Parameters &trackParameters) const
{
    const pandora::CartesianVector &referencePoint(pHelix->GetReferencePoint());

    // First project to endcap
    float minGenericTime(std::numeric_limits<float>::max());
    bool isProjectedToEndCap(true);

    pandora::CartesianVector bestECalProjection(0.f, 0.f, 0.f);
    (void) pHelix->GetPointInZ(static_cast<float>(signPz) * m_eCalEndCapInnerZ, referencePoint, bestECalProjection, minGenericTime);

    // Then project to barrel surface(s)
    pandora::CartesianVector barrelProjection(0.f, 0.f, 0.f);

    if (m_eCalBarrelInnerSymmetry > 0)
    {
        // Polygon
        float twopi_n = 2. * M_PI / (static_cast<float>(m_eCalBarrelInnerSymmetry));

        for (int i = 0; i < m_eCalBarrelInnerSymmetry; ++i)
        {
            float genericTime(std::numeric_limits<float>::max());
            const float phi(twopi_n * static_cast<float>(i) + m_eCalBarrelInnerPhi0);

            const pandora::StatusCode statusCode(pHelix->GetPointInXY(m_eCalBarrelInnerR * std::cos(phi), m_eCalBarrelInnerR * std::sin(phi),
                std::cos(phi + 0.5 * M_PI), std::sin(phi + 0.5 * M_PI), referencePoint, barrelProjection, genericTime));

            if ((pandora::STATUS_CODE_SUCCESS == statusCode) && (genericTime < minGenericTime))
            {
                minGenericTime = genericTime;
                isProjectedToEndCap = false;
                bestECalProjection = barrelProjection;
            }
        }
    }
    else
    {
        // Cylinder
        float genericTime(std::numeric_limits<float>::max());
        const pandora::StatusCode statusCode(pHelix->GetPointOnCircle(m_eCalBarrelInnerR, referencePoint, barrelProjection, genericTime));

        if ((pandora::STATUS_CODE_SUCCESS == statusCode) && (genericTime < minGenericTime))
        {
            minGenericTime = genericTime;
            isProjectedToEndCap = false;
            bestECalProjection = barrelProjection;
        }
    }

    if (pandora::CartesianVector(0.f, 0.f, 0.f) == bestECalProjection)
        throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_INITIALIZED);

    trackParameters.m_trackStateAtCalorimeter = pandora::TrackState(bestECalProjection, pHelix->GetExtrapolatedMomentum(bestECalProjection));
    trackParameters.m_isProjectedToEndCap = isProjectedToEndCap;

    // Convert generic time (length from reference point to intersection, divided by momentum) into nanoseconds
    const float particleMass(trackParameters.m_mass.Get());
    const float particleEnergy(std::sqrt(particleMass * particleMass + trackParameters.m_momentumAtDca.Get().GetMagnitudeSquared()));
    trackParameters.m_timeAtCalorimeter = minGenericTime * particleEnergy / 300.f;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackCreator::TrackReachesECAL(const EVENT::Track *const pTrack, PandoraApi::Track::Parameters &trackParameters) const
{
    // Calculate hit position information
    float hitZMin(std::numeric_limits<float>::max());
    float hitZMax(-std::numeric_limits<float>::max());
    float hitOuterR(-std::numeric_limits<float>::max());

    int nTpcHits(0);
    int nFtdHits(0);
    int maxOccupiedFtdLayer(0);

    const EVENT::TrackerHitVec &trackerHitVec(pTrack->getTrackerHits());
    const unsigned int nTrackHits(trackerHitVec.size());

    for (unsigned int i = 0; i < nTrackHits; ++i)
    {
        const float x(static_cast<float>(trackerHitVec[i]->getPosition()[0]));
        const float y(static_cast<float>(trackerHitVec[i]->getPosition()[1]));
        const float z(static_cast<float>(trackerHitVec[i]->getPosition()[2]));
        const float r(std::sqrt(x * x + y * y));

        if (z > hitZMax)
            hitZMax = z;

        if (z < hitZMin)
            hitZMin = z;

        if (r > hitOuterR)
            hitOuterR = r;

        if ((r > m_tpcInnerR) && (r < m_tpcOuterR) && (std::fabs(z) <= m_tpcZmax))
        {
            nTpcHits++;
            continue;
        }

        for (unsigned int j = 0; j < m_nFtdLayers; ++j)
        {
            if ((r > m_ftdInnerRadii[j]) && (r < m_ftdOuterRadii[j]) &&
                (std::fabs(z) - m_settings.m_reachesECalFtdZMaxDistance < m_ftdZPositions[j]) &&
                (std::fabs(z) + m_settings.m_reachesECalFtdZMaxDistance > m_ftdZPositions[j]))
            {
                if (static_cast<int>(j) > maxOccupiedFtdLayer)
                    maxOccupiedFtdLayer = static_cast<int>(j);

                nFtdHits++;
                break;
            }
        }
    }

    // Look to see if there are hits in etd or set, implying track has reached edge of ecal
    if ((hitOuterR > m_minSetRadius) || (hitZMax > m_minEtdZPosition))
    {
        trackParameters.m_reachesCalorimeter = true;
        return;
    }

    // Require sufficient hits in tpc or ftd, then compare extremal hit positions with tracker dimensions
    if ((nTpcHits >= m_settings.m_reachesECalNTpcHits) || (nFtdHits >= m_settings.m_reachesECalNFtdHits))
    {
        if ((hitOuterR - m_tpcOuterR > m_settings.m_reachesECalTpcOuterDistance) ||
            (std::fabs(hitZMax) - m_tpcZmax > m_settings.m_reachesECalTpcZMaxDistance) ||
            (std::fabs(hitZMin) - m_tpcZmax > m_settings.m_reachesECalTpcZMaxDistance) ||
            (maxOccupiedFtdLayer >= m_settings.m_reachesECalMinFtdLayer))
        {
            trackParameters.m_reachesCalorimeter = true;
            return;
        }
    }

    // If track is lowpt, it may curl up and end inside tpc inner radius
    const pandora::CartesianVector &momentumAtDca(trackParameters.m_momentumAtDca.Get());
    const float cosAngleAtDca(std::fabs(momentumAtDca.GetZ()) / momentumAtDca.GetMagnitude());
    const float pX(momentumAtDca.GetX()), pY(momentumAtDca.GetY());
    const float pT(std::sqrt(pX * pX + pY * pY));

    if ((cosAngleAtDca > m_cosTpc) || (pT < m_settings.m_curvatureToMomentumFactor * m_bField * m_tpcOuterR))
    {
        trackParameters.m_reachesCalorimeter = true;
        return;
    }

    trackParameters.m_reachesCalorimeter = false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackCreator::DefineTrackPfoUsage(const EVENT::Track *const pTrack, PandoraApi::Track::Parameters &trackParameters) const
{
    bool canFormPfo(false);
    bool canFormClusterlessPfo(false);

    if (trackParameters.m_reachesCalorimeter.Get() && !this->IsParent(pTrack))
    {
        const float d0(std::fabs(pTrack->getD0())), z0(std::fabs(pTrack->getZ0()));

        EVENT::TrackerHitVec trackerHitvec(pTrack->getTrackerHits());
        float rInner(std::numeric_limits<float>::max()), zMin(std::numeric_limits<float>::max());

        for (EVENT::TrackerHitVec::const_iterator iter = trackerHitvec.begin(), iterEnd = trackerHitvec.end(); iter != iterEnd; ++iter)
        {
            const double *pPosition((*iter)->getPosition());
            const float x(pPosition[0]), y(pPosition[1]), absoluteZ(std::fabs(pPosition[2]));
            const float r(std::sqrt(x * x + y * y));

            if (r < rInner)
                rInner = r;

            if (absoluteZ < zMin)
                zMin = absoluteZ;
        }

        if (this->PassesQualityCuts(pTrack, trackParameters))
        {
            const pandora::CartesianVector &momentumAtDca(trackParameters.m_momentumAtDca.Get());
            const float pX(momentumAtDca.GetX()), pY(momentumAtDca.GetY()), pZ(momentumAtDca.GetZ());
            const float pT(std::sqrt(pX * pX + pY * pY));

            const float zCutForNonVertexTracks(m_tpcInnerR * std::fabs(pZ / pT) + m_settings.m_zCutForNonVertexTracks);
            const bool passRzQualityCuts((zMin < zCutForNonVertexTracks) && (rInner < m_tpcInnerR + m_settings.m_maxTpcInnerRDistance));

            const bool isV0(this->IsV0(pTrack));
            const bool isDaughter(this->IsDaughter(pTrack));

            // Decide whether track can be associated with a pandora cluster and used to form a charged PFO
            if ((d0 < m_settings.m_d0TrackCut) && (z0 < m_settings.m_z0TrackCut) && (rInner < m_tpcInnerR + m_settings.m_maxTpcInnerRDistance))
            {
                canFormPfo = true;
            }
            else if (passRzQualityCuts && (0 != m_settings.m_usingNonVertexTracks))
            {
                canFormPfo = true;
            }
            else if (isV0 || isDaughter)
            {
                canFormPfo = true;
            }

            // Decide whether track can be used to form a charged PFO, even if track fails to be associated with a pandora cluster
            const float particleMass(trackParameters.m_mass.Get());
            const float trackEnergy(std::sqrt(momentumAtDca.GetMagnitudeSquared() + particleMass * particleMass));

            if ((0 != m_settings.m_usingUnmatchedVertexTracks) && (trackEnergy < m_settings.m_unmatchedVertexTrackMaxEnergy))
            {
                if ((d0 < m_settings.m_d0UnmatchedVertexTrackCut) && (z0 < m_settings.m_z0UnmatchedVertexTrackCut) &&
                    (rInner < m_tpcInnerR + m_settings.m_maxTpcInnerRDistance))
                {
                    canFormClusterlessPfo = true;
                }
                else if (passRzQualityCuts && (0 != m_settings.m_usingNonVertexTracks) && (0 != m_settings.m_usingUnmatchedNonVertexTracks))
                {
                    canFormClusterlessPfo = true;
                }
                else if (isV0 || isDaughter)
                {
                    canFormClusterlessPfo = true;
                }
            }
        }
        else if (this->IsDaughter(pTrack) || this->IsV0(pTrack))
        {
            streamlog_out(WARNING) << "Recovering daughter or v0 track " << trackParameters.m_momentumAtDca.Get().GetMagnitude() << std::endl;
            canFormPfo = true;
        }
    }

    trackParameters.m_canFormPfo = canFormPfo;
    trackParameters.m_canFormClusterlessPfo = canFormClusterlessPfo;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TrackCreator::PassesQualityCuts(const EVENT::Track *const pTrack, const PandoraApi::Track::Parameters &trackParameters) const
{
    // First simple sanity checks
    if (trackParameters.m_trackStateAtCalorimeter.Get().GetPosition().GetMagnitude() < m_settings.m_minTrackECalDistanceFromIp)
        return false;

    if (pTrack->getOmega() == 0.f)
    {
        streamlog_out(ERROR) << "Track has Omega = 0 " << std::endl;
        return false;
    }

    // Check momentum uncertainty is reasonable to use track
    const pandora::CartesianVector &momentumAtDca(trackParameters.m_momentumAtDca.Get());
    const float sigmaPOverP(std::sqrt(pTrack->getCovMatrix()[5]) / std::fabs(pTrack->getOmega()));

    if (sigmaPOverP > m_settings.m_maxTrackSigmaPOverP)
    {
        streamlog_out(WARNING) << " Dropping track : " << momentumAtDca.GetMagnitude() << "+-" << sigmaPOverP * (momentumAtDca.GetMagnitude())
                               << " chi2 = " <<  pTrack->getChi2() << " " << pTrack->getNdf()
                               << " from " << pTrack->getTrackerHits().size() << std::endl;
        return false;
    }

    // Require reasonable number of TPC hits 
    if (momentumAtDca.GetMagnitude() > m_settings.m_minMomentumForTrackHitChecks)
    {
        const float pX(fabs(momentumAtDca.GetX()));
        const float pY(fabs(momentumAtDca.GetY()));
        const float pZ(fabs(momentumAtDca.GetZ()));
        const float pT(std::sqrt(pX * pX + pY * pY));
        const float rInnermostHit(pTrack->getRadiusOfInnermostHit());

        if ((0.f == pT) || (0.f == pZ) || (rInnermostHit == m_tpcOuterR))
        {
            streamlog_out(ERROR) << "Invalid track parameter, pT " << pT << ", pZ " << pZ << ", rInnermostHit " << rInnermostHit << std::endl;
            return false;
        }

        float nExpectedTpcHits(0.);

        if (pZ < m_tpcZmax / m_tpcOuterR * pT)
        {
            const float innerExpectedHitRadius(std::max(m_tpcInnerR, rInnermostHit));
            const float frac((m_tpcOuterR - innerExpectedHitRadius) / (m_tpcOuterR - m_tpcInnerR));
            nExpectedTpcHits = m_tpcMaxRow * frac;
        }

        if ((pZ <= m_tpcZmax / m_tpcInnerR * pT) && (pZ >= m_tpcZmax / m_tpcOuterR * pT))
        {
            const float innerExpectedHitRadius(std::max(m_tpcInnerR, rInnermostHit));
            const float frac((m_tpcZmax * pT / pZ - innerExpectedHitRadius) / (m_tpcOuterR - innerExpectedHitRadius));
            nExpectedTpcHits = frac * m_tpcMaxRow;
        }

        // TODO Get TPC membrane information from GEAR when available
        if (std::fabs(pZ) / momentumAtDca.GetMagnitude() < m_settings.m_tpcMembraneMaxZ / m_tpcInnerR)
            nExpectedTpcHits = 0;

        const EVENT::IntVec &hitsBySubdetector(pTrack->getSubdetectorHitNumbers());

        //fg: hit numbers are now given in different order wrt LOI:  
        // trk->subdetectorHitNumbers()[ 2 * ILDDetID::TPC - 1 ] =  hitsInFit ;  
        // trk->subdetectorHitNumbers()[ 2 * ILDDetID::TPC - 2 ] =  hitCount ;  
        // ---- use hitsInFit :
        const int nTpcHits = hitsBySubdetector[ 2 * lcio::ILDDetID::TPC - 1 ];
        const int nFtdHits = hitsBySubdetector[ 2 * lcio::ILDDetID::FTD - 1 ];

        const int minTpcHits = static_cast<int>(nExpectedTpcHits * m_settings.m_minTpcHitFractionOfExpected);

        if ((nTpcHits < minTpcHits) && (nFtdHits < m_settings.m_minFtdHitsForTpcHitFraction))
        {
            streamlog_out(WARNING) << " Dropping track : " << momentumAtDca.GetMagnitude() << " Number of TPC hits = " << nTpcHits
                                   << " < " << minTpcHits << " nftd = " << nFtdHits  << std::endl;
            return false;
        }
    }

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

TrackCreator::Settings::Settings() :
    m_shouldFormTrackRelationships(1),
    m_minTrackHits(5),
    m_minFtdTrackHits(0),
    m_maxTrackHits(5000.f),
    m_useOldTrackStateCalculation(0),
    m_d0TrackCut(50.f),
    m_z0TrackCut(50.f),
    m_usingNonVertexTracks(1),
    m_usingUnmatchedNonVertexTracks(0),
    m_usingUnmatchedVertexTracks(1),
    m_unmatchedVertexTrackMaxEnergy(5.f),
    m_d0UnmatchedVertexTrackCut(5.f),
    m_z0UnmatchedVertexTrackCut(5.f),
    m_zCutForNonVertexTracks(250.f),
    m_reachesECalNTpcHits(11),
    m_reachesECalNFtdHits(4),
    m_reachesECalTpcOuterDistance(-100.f),
    m_reachesECalMinFtdLayer(9),
    m_reachesECalTpcZMaxDistance(-50.f),
    m_reachesECalFtdZMaxDistance(1.f),
    m_curvatureToMomentumFactor(0.3f / 2000.f),
    m_minTrackECalDistanceFromIp(100.f),
    m_maxTrackSigmaPOverP(0.15f),
    m_minMomentumForTrackHitChecks(1.f),
    m_tpcMembraneMaxZ(10.f),
    m_maxTpcInnerRDistance(50.f),
    m_minTpcHitFractionOfExpected(0.2f),
    m_minFtdHitsForTpcHitFraction(2)
{
}
