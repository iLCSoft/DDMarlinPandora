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

#include "DDTrackCreatorBase.h"
#include "Pandora/PdgTable.h"
#include <LCObjects/LCTrack.h>

#include <MarlinTrk/Factory.h>
#include <MarlinTrk/IMarlinTrack.h>

#include "DD4hep/Detector.h"
#include "DD4hep/DD4hepUnits.h"
#include "DDRec/DetectorData.h"

#include <algorithm>
#include <cmath>
#include <limits>

//forward declaration
std::vector<double> getTrackingRegionExtent();

DDTrackCreatorBase::DDTrackCreatorBase(const Settings &settings, const pandora::Pandora *const pPandora) :
    m_settings(settings),
    m_pandora(*pPandora),
    m_trackVector(0),
    m_v0TrackList( TrackList() ),
    m_parentTrackList( TrackList() ),
    m_daughterTrackList( TrackList() ),
    m_trackToPidMap( TrackToPidMap() ),
    m_minimalTrackStateRadiusSquared( 0.f )
{

    const float ecalInnerR = settings.m_eCalBarrelInnerR;
    const float tsTolerance = settings.m_trackStateTolerance;
    m_minimalTrackStateRadiusSquared = (ecalInnerR-tsTolerance)*(ecalInnerR-tsTolerance);
    //wrap in shared_ptr with a dummy destructor
    m_trackingSystem =
      std::shared_ptr<MarlinTrk::IMarlinTrkSystem>( MarlinTrk::Factory::createMarlinTrkSystem(settings.m_trackingSystemName,
                                                                                              nullptr, ""),
                                                    [](MarlinTrk::IMarlinTrkSystem*){} );
    m_trackingSystem->init();
    m_encoder = std::make_shared<UTIL::BitField64>( lcio::LCTrackerCellID::encoding_string() );
    m_lcTrackFactory = std::make_shared<lc_content::LCTrackFactory>();
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
                                PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::SetTrackParentDaughterRelationship(m_pandora,
                                    pTrack, trackVec[jTrack]));
                            }
                        }

                        // Make track sibling relationships
                        else
                        {
                            for (unsigned int jTrack = iTrack + 1; jTrack < nTracks; ++jTrack)
                            {
                                PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::SetTrackSiblingRelationship(m_pandora,
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
                                PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::SetTrackParentDaughterRelationship(m_pandora,
                                    pTrack, trackVec[jTrack]));
                            }
                        }

                        // Make track sibling relationships
                        else
                        {
                            for (unsigned int jTrack = iTrack + 1; jTrack < nTracks; ++jTrack)
                            {
                                PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::SetTrackSiblingRelationship(m_pandora,
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
                            PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::SetTrackSiblingRelationship(m_pandora,
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

    // Convert generic time (length from reference point to intersection, divided by momentum) into nanoseconds
    const float minGenericTime(this->CalculateTrackTimeAtCalorimeter(pTrack));
    const float particleMass(trackParameters.m_mass.Get());
    const float particleEnergy(std::sqrt(particleMass * particleMass + trackParameters.m_momentumAtDca.Get().GetMagnitudeSquared()));
    trackParameters.m_timeAtCalorimeter = minGenericTime * particleEnergy / 299.792f;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DDTrackCreatorBase::GetTrackStatesAtCalo( EVENT::Track *track,
                                               lc_content::LCTrackParameters& trackParameters ){

  if( not trackParameters.m_reachesCalorimeter.Get() ) {
      streamlog_out(DEBUG5) << "Track does not reach the ECal" <<std::endl;
    return;
  }

  const TrackState *trackAtCalo = track->getTrackState(TrackState::AtCalorimeter);
  if( not trackAtCalo ) {
      streamlog_out(DEBUG5) << "Track does not have a trackState at calorimeter" <<std::endl;
      streamlog_out(DEBUG3) << toString(track) << std::endl;
      return;
  }

  streamlog_out(DEBUG3) << "Original" << toString(trackAtCalo) << std::endl;

  const auto* tsPosition = trackAtCalo->getReferencePoint();

  if( std::fabs(tsPosition[2]) <  getTrackingRegionExtent()[2] ) {
      streamlog_out(DEBUG5) << "Original trackState is at Barrel" << std::endl;
      pandora::InputTrackState pandoraTrackState;
      this->CopyTrackState( trackAtCalo, pandoraTrackState );
      trackParameters.m_trackStates.push_back( pandoraTrackState );
  } else { // if track state is in endcap we do not repeat track state calculation, because the barrel cannot be hit
      streamlog_out(DEBUG5) << "Original track state is at Endcap" << std::endl;
      pandora::InputTrackState pandoraTrackState;
      this->CopyTrackState( trackAtCalo, pandoraTrackState );
      trackParameters.m_trackStates.push_back( pandoraTrackState );
      return;
  }

  auto marlintrk = std::unique_ptr<MarlinTrk::IMarlinTrack>(m_trackingSystem->createTrack());
  const EVENT::TrackerHitVec& trkHits = track->getTrackerHits();
  const int nHitsTrack = trkHits.size();

  for (int iHit = 0; iHit < nHitsTrack; ++iHit) {
    EVENT::TrackerHit* trkHit = trkHits[iHit] ;
    if( UTIL::BitSet32( trkHit->getType() )[ UTIL::ILDTrkHitTypeBit::COMPOSITE_SPACEPOINT ]   ){ //it is a composite spacepoint
      //Split it up and add both hits to the MarlinTrk
      const EVENT::LCObjectVec& rawObjects = trkHit->getRawHits();
      for( unsigned k=0; k< rawObjects.size(); k++ ){
	EVENT::TrackerHit* rawHit = static_cast< EVENT::TrackerHit* >( rawObjects[k] );
	if( marlintrk->addHit( rawHit ) != MarlinTrk::IMarlinTrack::success ){
	  streamlog_out(DEBUG4) << "DDTrackCreatorBase::GetTrackStatesAtCalo failed to add strip hit " << *rawHit << std::endl;
	}
      }
    } else {
      if( marlintrk->addHit(trkHits[iHit])  != MarlinTrk::IMarlinTrack::success  )
	streamlog_out(DEBUG4) << "DDTrackCreatorBase::GetTrackStatesAtCalo failed to add tracker hit " << *trkHit<< std::endl;
    }
  }

  bool tanL_is_positive = trackAtCalo->getTanLambda()>0;

  auto trackState = TrackStateImpl(*trackAtCalo);

  int return_error  = marlintrk->initialise(trackState, m_settings.m_bField, MarlinTrk::IMarlinTrack::modeForward);
  if (return_error != MarlinTrk::IMarlinTrack::success ) {
    streamlog_out(DEBUG4) << "DDTrackCreatorBase::GetTrackStatesAtCalo failed to initialize track for endcap track : " << std::endl ;
    return ;
  }

  double chi2 = -DBL_MAX;
  int ndf = 0;

  TrackStateImpl trackStateAtCaloEndcap;

  unsigned ecal_endcap_face_ID = lcio::ILDDetID::ECAL_ENDCAP;
  int detElementID = 0;
  m_encoder->reset();  // reset to 0
  (*m_encoder)[lcio::LCTrackerCellID::subdet()] = ecal_endcap_face_ID;
  (*m_encoder)[lcio::LCTrackerCellID::side()] = tanL_is_positive ? lcio::ILDDetID::fwd : lcio::ILDDetID::bwd;
  (*m_encoder)[lcio::LCTrackerCellID::layer()]  = 0;

  return_error = marlintrk->propagateToLayer(m_encoder->lowWord(), trackStateAtCaloEndcap, chi2, ndf,
                                             detElementID, MarlinTrk::IMarlinTrack::modeForward );
  streamlog_out(DEBUG5) << "Found trackState at endcap? Error code: " << return_error  << std::endl;

  if (return_error == MarlinTrk::IMarlinTrack::success ) {
      streamlog_out(DEBUG3) << "Endcap" << toString(&trackStateAtCaloEndcap) << std::endl;
      const auto* tsEP = trackStateAtCaloEndcap.getReferencePoint();
      const double radSquared = ( tsEP[0]*tsEP[0] + tsEP[1]*tsEP[1] );
      if( radSquared < m_minimalTrackStateRadiusSquared ) {
          streamlog_out(DEBUG5) << "new track state is below tolerance radius" << std::endl;
          return;
      }
      //for curling tracks the propagated track has the wrong z0 whereas it should be 0. really
      if( std::abs( trackStateAtCaloEndcap.getZ0() ) >
          std::abs( 2.*M_PI/trackStateAtCaloEndcap.getOmega() * trackStateAtCaloEndcap.getTanLambda() ) ){
          trackStateAtCaloEndcap.setZ0( 0. );
      }
      streamlog_out(DEBUG5) << "new track state at endcap accepted" << std::endl;
      pandora::InputTrackState pandoraAtEndcap;
      this->CopyTrackState( &trackStateAtCaloEndcap, pandoraAtEndcap );
      trackParameters.m_trackStates.push_back( pandoraAtEndcap );
  }

  return;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float DDTrackCreatorBase::CalculateTrackTimeAtCalorimeter(const EVENT::Track *const pTrack) const
{
    const pandora::Helix helix(pTrack->getPhi(), pTrack->getD0(), pTrack->getZ0(), pTrack->getOmega(), pTrack->getTanLambda(), m_settings.m_bField);
    const pandora::CartesianVector &referencePoint(helix.GetReferencePoint());

    // First project to endcap
    float minGenericTime(std::numeric_limits<float>::max());

    pandora::CartesianVector bestECalProjection(0.f, 0.f, 0.f);
    const int signPz((helix.GetMomentum().GetZ() > 0.f) ? 1 : -1);
    (void) helix.GetPointInZ(static_cast<float>(signPz) * m_settings.m_eCalEndCapInnerZ, referencePoint, bestECalProjection, minGenericTime);

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

            const pandora::StatusCode statusCode(helix.GetPointInXY(m_settings.m_eCalBarrelInnerR * std::cos(phi), m_settings.m_eCalBarrelInnerR * std::sin(phi),
                std::cos(phi + 0.5 * M_PI), std::sin(phi + 0.5 * M_PI), referencePoint, barrelProjection, genericTime));

            if ((pandora::STATUS_CODE_SUCCESS == statusCode) && (genericTime < minGenericTime))
            {
                minGenericTime = genericTime;
                bestECalProjection = barrelProjection;
            }
        }
    }
    else
    {
        // Cylinder
        float genericTime(std::numeric_limits<float>::max());
        const pandora::StatusCode statusCode(helix.GetPointOnCircle(m_settings.m_eCalBarrelInnerR, referencePoint, barrelProjection, genericTime));

        if ((pandora::STATUS_CODE_SUCCESS == statusCode) && (genericTime < minGenericTime))
        {
            minGenericTime = genericTime;
            bestECalProjection = barrelProjection;
        }
    }

    if (bestECalProjection.GetMagnitudeSquared() < std::numeric_limits<float>::epsilon())
        throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_INITIALIZED);

    return minGenericTime;
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

DDTrackCreatorBase::Settings::Settings() :
    m_trackCollections( StringVector() ),
    m_kinkVertexCollections( StringVector() ),
    m_prongVertexCollections( StringVector() ),
    m_splitVertexCollections( StringVector() ),
    m_v0VertexCollections( StringVector() ),
    m_prongSplitVertexCollections(StringVector()),
    m_shouldFormTrackRelationships(1),
    m_minTrackHits(5),
    m_minFtdTrackHits(0),
    m_maxTrackHits(5000.f),
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
    m_minFtdHitsForBarrelTrackerHitFraction(2),
    m_trackStateTolerance(0.f),
    m_trackingSystemName("DDKalTest"),
    m_bField(0.f),
    m_eCalBarrelInnerSymmetry(0),
    m_eCalBarrelInnerPhi0(0.f),
    m_eCalBarrelInnerR(0.f),
    m_eCalEndCapInnerZ(0.f)
{
}
