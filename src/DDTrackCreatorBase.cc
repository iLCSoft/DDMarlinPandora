/**
 *  @file   DDMarlinPandora/src/DDTrackCreatorBase.cc
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
#include "UTIL/Operators.h"

#include "DDPandoraPFANewProcessor.h"
#include "DDTrackCreatorBase.h"
#include "Pandora/PdgTable.h"

#include <algorithm>
#include <cmath>
#include <limits>

#include "DD4hep/LCDD.h"
#include "DD4hep/DD4hepUnits.h"
#include "DDRec/DetectorData.h"



//forward declarations. See in DDPandoraPFANewProcessor.cc
double getFieldFromLCDD(); 
DD4hep::DDRec::LayeredCalorimeterData * getExtension(std::string detectorName);

std::vector<double> getTrackingRegionExtent();

DDTrackCreatorBase::DDTrackCreatorBase(const Settings &settings, const pandora::Pandora *const pPandora) :
    m_settings(settings),
    m_pPandora(pPandora)
{
  
}

//------------------------------------------------------------------------------------------------------------------------------------------

DDTrackCreatorBase::~DDTrackCreatorBase()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode DDTrackCreatorBase::CreateTrackAssociations(const EVENT::LCEvent *const pLCEvent)
{
    PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, this->ExtractKinks(pLCEvent));
    PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, this->ExtractProngsAndSplits(pLCEvent));
    PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, this->ExtractV0s(pLCEvent));

    return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode DDTrackCreatorBase::ExtractKinks(const EVENT::LCEvent *const pLCEvent)
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

pandora::StatusCode DDTrackCreatorBase::ExtractProngsAndSplits(const EVENT::LCEvent *const pLCEvent)
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

pandora::StatusCode DDTrackCreatorBase::ExtractV0s(const EVENT::LCEvent *const pLCEvent)
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

bool DDTrackCreatorBase::IsConflictingRelationship(const EVENT::TrackVec &trackVec) const
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

void DDTrackCreatorBase::GetTrackStates(const EVENT::Track *const pTrack, PandoraApi::Track::Parameters &trackParameters) const
{
    if (m_settings.m_useOldTrackStateCalculation > 0)
        return this->GetTrackStatesOld(pTrack, trackParameters);

    const TrackState *pTrackState = pTrack->getTrackState(TrackState::AtIP);

    if (!pTrackState)
        throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_INITIALIZED);

    const double pt(m_settings.m_bField * 2.99792e-4 / std::fabs(pTrackState->getOmega()));
    trackParameters.m_momentumAtDca = pandora::CartesianVector(std::cos(pTrackState->getPhi()), std::sin(pTrackState->getPhi()), pTrackState->getTanLambda()) * pt;

    this->CopyTrackState(pTrack->getTrackState(TrackState::AtFirstHit), trackParameters.m_trackStateAtStart);

    //fg: curling TPC tracks have pointers to track segments stored -> need to get track states from last segment!
    const EVENT::Track *pEndTrack = (pTrack->getTracks().empty() ?  pTrack  :  pTrack->getTracks().back());

    this->CopyTrackState(pEndTrack->getTrackState(TrackState::AtLastHit), trackParameters.m_trackStateAtEnd);
    this->CopyTrackState(pEndTrack->getTrackState(TrackState::AtCalorimeter), trackParameters.m_trackStateAtCalorimeter);

    trackParameters.m_isProjectedToEndCap = ((std::fabs(trackParameters.m_trackStateAtCalorimeter.Get().GetPosition().GetZ()) < m_settings.m_eCalEndCapInnerZ) ? false : true);

    // TODO minGenericTime * particleEnergy / 300.f;
    trackParameters.m_timeAtCalorimeter = 0.f;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DDTrackCreatorBase::CopyTrackState(const TrackState *const pTrackState, pandora::InputTrackState &inputTrackState) const
{
    if (!pTrackState)
        throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_INITIALIZED);

    const double pt(m_settings.m_bField * 2.99792e-4 / std::fabs(pTrackState->getOmega()));

    const double px(pt * std::cos(pTrackState->getPhi()));
    const double py(pt * std::sin(pTrackState->getPhi()));
    const double pz(pt * pTrackState->getTanLambda());

    const double xs(pTrackState->getReferencePoint()[0]);
    const double ys(pTrackState->getReferencePoint()[1]);
    const double zs(pTrackState->getReferencePoint()[2]);

    inputTrackState = pandora::TrackState(xs, ys, zs, px, py, pz);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DDTrackCreatorBase::GetTrackStatesOld(const EVENT::Track *const pTrack, PandoraApi::Track::Parameters &trackParameters) const
{
    pandora::Helix *pHelixFit = new pandora::Helix(pTrack->getPhi(), pTrack->getD0(), pTrack->getZ0(), pTrack->getOmega(), pTrack->getTanLambda(), m_settings.m_bField);
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

void DDTrackCreatorBase::GetECalProjectionOld(const pandora::Helix *const pHelix, const int signPz, PandoraApi::Track::Parameters &trackParameters) const
{
    const pandora::CartesianVector &referencePoint(pHelix->GetReferencePoint());

    // First project to endcap
    float minGenericTime(std::numeric_limits<float>::max());
    bool isProjectedToEndCap(true);

    pandora::CartesianVector bestECalProjection(0.f, 0.f, 0.f);
    (void) pHelix->GetPointInZ(static_cast<float>(signPz) * m_settings.m_eCalEndCapInnerZ, referencePoint, bestECalProjection, minGenericTime);

    // Then project to barrel surface(s)
    pandora::CartesianVector barrelProjection(0.f, 0.f, 0.f);

    if (m_settings.m_eCalBarrelInnerSymmetry > 0)
    {
        // Polygon
        float twopi_n = 2. * M_PI / (static_cast<float>(m_settings.m_eCalBarrelInnerSymmetry));

        for (int i = 0; i < m_settings.m_eCalBarrelInnerSymmetry; ++i)
        {
            float genericTime(std::numeric_limits<float>::max());
            const float phi(twopi_n * static_cast<float>(i) + m_settings.m_eCalBarrelInnerPhi0);

            const pandora::StatusCode statusCode(pHelix->GetPointInXY(m_settings.m_eCalBarrelInnerR * std::cos(phi), m_settings.m_eCalBarrelInnerR * std::sin(phi),
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
        const pandora::StatusCode statusCode(pHelix->GetPointOnCircle(m_settings.m_eCalBarrelInnerR, referencePoint, barrelProjection, genericTime));

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


DDTrackCreatorBase::Settings::Settings() :
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
    m_reachesECalNBarrelTrackerHits(11),
    m_reachesECalNFtdHits(4),
    m_reachesECalBarrelTrackerOuterDistance(-100.f),
    m_reachesECalMinFtdLayer(9),
    m_reachesECalBarrelTrackerZMaxDistance(-50.f),
    m_reachesECalFtdZMaxDistance(1.f),
    m_curvatureToMomentumFactor(0.3f / 2000.f),
    m_minTrackECalDistanceFromIp(100.f),
    m_maxTrackSigmaPOverP(0.15f),
    m_minMomentumForTrackHitChecks(1.f),
    m_maxBarrelTrackerInnerRDistance(50.f),
    m_minBarrelTrackerHitFractionOfExpected(0.2f),
    m_minFtdHitsForBarrelTrackerHitFraction(2)
{
}
