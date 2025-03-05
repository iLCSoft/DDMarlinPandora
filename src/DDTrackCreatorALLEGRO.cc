/**
 *  @file   DDMarlinPandora/src/DDTrackCreatorALLEGRO.cc
 * 
 *  @brief  Implementation of the track creator class for a ALLEGRO all silicon tracker.
 * 
 *  $Log: $
 */

#include "DDTrackCreatorALLEGRO.h"

#include "marlin/Global.h"
#include "marlin/Processor.h"

#include "EVENT/LCCollection.h"
#include "EVENT/ReconstructedParticle.h"
#include "EVENT/Vertex.h"
#include "UTIL/ILDConf.h"
#include "UTIL/Operators.h"

#include "Pandora/PdgTable.h"
#include "LCObjects/LCTrack.h"

#include "DD4hep/Detector.h"
#include "DD4hep/DD4hepUnits.h"
#include "DDRec/DetectorData.h"
#include "DD4hep/DetType.h"
#include "DD4hep/DetectorSelector.h"

#include <algorithm>
#include <cmath>
#include <limits>

//forward declarations. See in DDPandoraPFANewProcessor.cc
std::vector<double> getTrackingRegionExtent();

DDTrackCreatorALLEGRO::DDTrackCreatorALLEGRO(const Settings &settings, const pandora::Pandora *const pPandora)
  : DDTrackCreatorBase(settings, pPandora),
    m_trackerInnerR( 0.f ),
    m_trackerOuterR( 0.f ),
    m_trackerZmax( 0.f ),
    m_cosTracker( 0.f ),
    m_endcapDiskInnerRadii( DoubleVector() ),
    m_endcapDiskOuterRadii( DoubleVector() ),
    m_endcapDiskZPositions( DoubleVector() ),
    m_nEndcapDiskLayers( 0 ),
    m_barrelTrackerLayers( 0 ),
    m_tanLambdaEndcapDisk( 0.f )

{
    // FIXME! AD: currently ignoring the tracker parameters since we are working with truth tracks...
    /*
    m_trackerInnerR = getTrackingRegionExtent()[0];
    m_trackerOuterR = getTrackingRegionExtent()[1];
    m_trackerZmax = getTrackingRegionExtent()[2];

    ///FIXME: Probably need to be something relating to last disk inner radius
    m_cosTracker = m_trackerZmax / std::sqrt(m_trackerZmax * m_trackerZmax + m_trackerInnerR * m_trackerInnerR);

    dd4hep::Detector & mainDetector = dd4hep::Detector::getInstance();

    //Maybe we need to veto the vertex? That was done in the ILD case
    const std::vector< dd4hep::DetElement>& barrelDets = dd4hep::DetectorSelector(mainDetector).detectors(  ( dd4hep::DetType::TRACKER | dd4hep::DetType::BARREL )) ;

    m_barrelTrackerLayers = 0;
    for (std::vector< dd4hep::DetElement>::const_iterator iter = barrelDets.begin(), iterEnd = barrelDets.end();iter != iterEnd; ++iter){
        try
        {
            dd4hep::rec::ZPlanarData * theExtension = 0;

            const dd4hep::DetElement& theDetector = *iter;
            theExtension = theDetector.extension<dd4hep::rec::ZPlanarData>();

            unsigned int N = theExtension->layers.size();
            m_barrelTrackerLayers=m_barrelTrackerLayers+N;

            streamlog_out( DEBUG2 ) << " Adding layers for barrel tracker from DD4hep for "<< theDetector.name()<< "- n layers: " << N<< " sum up to now: "<<m_barrelTrackerLayers<<std::endl;
        } catch (std::runtime_error &exception){
            streamlog_out(WARNING) << "DDTrackCreatorALLEGRO exception during Barrel Tracker layer sum for "<<const_cast<dd4hep::DetElement&>(*iter).name()<<" : " << exception.what() << std::endl;
        }
    }

    m_nEndcapDiskLayers=0;
    m_endcapDiskInnerRadii.clear();
    m_endcapDiskOuterRadii.clear();

    //Instead of gear, loop over a provided list of forward (read: endcap) tracking detectors. For ILD this would be FTD
     const std::vector< dd4hep::DetElement>& endcapDets = dd4hep::DetectorSelector(mainDetector).detectors(  ( dd4hep::DetType::TRACKER | dd4hep::DetType::ENDCAP )) ;

     for (std::vector< dd4hep::DetElement>::const_iterator iter = endcapDets.begin(), iterEnd = endcapDets.end();iter != iterEnd; ++iter){
        try
        {
            dd4hep::rec::ZDiskPetalsData * theExtension = 0;

            const dd4hep::DetElement& theDetector = *iter;
            theExtension = theDetector.extension<dd4hep::rec::ZDiskPetalsData>();

            unsigned int N = theExtension->layers.size();
            streamlog_out( DEBUG2 ) << " Filling FTD-like parameters from DD4hep for "<< theDetector.name()<< "- n layers: " << N<< std::endl;

            for(unsigned int i = 0; i < N; ++i)
            {
                dd4hep::rec::ZDiskPetalsData::LayerLayout thisLayer  = theExtension->layers[i];

                // Create a disk to represent even number petals front side
                //FIXME! VERIFY THAT TIS MAKES SENSE!
                m_endcapDiskInnerRadii.push_back(thisLayer.distanceSensitive/dd4hep::mm);
                m_endcapDiskOuterRadii.push_back(thisLayer.distanceSensitive/dd4hep::mm+thisLayer.lengthSensitive/dd4hep::mm);

                // Take the mean z position of the staggered petals
                const double zpos(thisLayer.zPosition/dd4hep::mm);
                m_endcapDiskZPositions.push_back(zpos);
                streamlog_out( DEBUG2 ) << "     layer " << i << " - mean z position = " << zpos << std::endl;
            }

            m_nEndcapDiskLayers = m_endcapDiskZPositions.size() ;
        } catch (std::runtime_error &exception){
            streamlog_out(WARNING) << "DDTrackCreatorALLEGRO exception during Forward Tracking Disk parameter construction for detector "<<const_cast<dd4hep::DetElement&>(*iter).name()<<" : " << exception.what() << std::endl;
        }
    }

    for (unsigned int iEndcapDiskLayer = 0; iEndcapDiskLayer < m_nEndcapDiskLayers; ++iEndcapDiskLayer)
    {
        if ((std::fabs(m_endcapDiskOuterRadii[iEndcapDiskLayer]) < std::numeric_limits<float>::epsilon()) ||
            (std::fabs(m_endcapDiskInnerRadii[iEndcapDiskLayer]) < std::numeric_limits<float>::epsilon()))
        {
            throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);
        }
    }

    m_tanLambdaEndcapDisk = m_endcapDiskZPositions[0] / m_endcapDiskOuterRadii[0];
    */
}

//------------------------------------------------------------------------------------------------------------------------------------------

DDTrackCreatorALLEGRO::~DDTrackCreatorALLEGRO()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode DDTrackCreatorALLEGRO::CreateTracks(EVENT::LCEvent *pLCEvent)
{
    for (StringVector::const_iterator iter = m_settings.m_trackCollections.begin(), iterEnd = m_settings.m_trackCollections.end();
        iter != iterEnd; ++iter)
    {
        try
        {
            const EVENT::LCCollection *pTrackCollection = pLCEvent->getCollection(*iter);
            for (int i = 0, iMax = pTrackCollection->getNumberOfElements(); i < iMax; ++i)
            {
                    EVENT::Track *pTrack = dynamic_cast<Track*>(pTrackCollection->getElementAt(i));

                    if (NULL == pTrack)
                        throw EVENT::Exception("Collection type mismatch");

		    streamlog_out(DEBUG0)<<" Warning! Ignoring expected number of hits and other hit number cuts. Should eventually change!"<<std::endl;

                    // Proceed to create the pandora track
                    lc_content::LCTrackParameters trackParameters;
                    trackParameters.m_d0 = pTrack->getD0();
                    trackParameters.m_z0 = pTrack->getZ0();
                    trackParameters.m_pParentAddress = pTrack;

                    // By default, assume tracks are charged pions
                    const float signedCurvature(pTrack->getOmega());
                    trackParameters.m_particleId = (signedCurvature > 0) ? pandora::PI_PLUS : pandora::PI_MINUS;
                    trackParameters.m_mass = pandora::PdgTable::GetParticleMass(pandora::PI_PLUS);

                    // Use particle id information from V0 and Kink finders
                    TrackToPidMap::const_iterator trackPIDiter = m_trackToPidMap.find(pTrack);

                    if(trackPIDiter != m_trackToPidMap.end())
                    {
                        trackParameters.m_particleId = trackPIDiter->second;
                        trackParameters.m_mass = pandora::PdgTable::GetParticleMass(trackPIDiter->second);
                    }

                    if (0.f != signedCurvature)
                        trackParameters.m_charge = static_cast<int>(signedCurvature / std::fabs(signedCurvature));

		    try
		      {
			this->GetTrackStates(pTrack, trackParameters);
			this->TrackReachesECAL(pTrack, trackParameters);
                        // FIXME! AD: this is disabled for ALLEGRO since I am not sure what to do here with truth tracks
			// GM: re-enable
			this->GetTrackStatesAtCalo(pTrack, trackParameters);
			this->DefineTrackPfoUsage(pTrack, trackParameters);

			PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::Track::Create(m_pandora, trackParameters, *m_lcTrackFactory));
			m_trackVector.push_back(pTrack);
		      }
		    catch (pandora::StatusCodeException &statusCodeException)
		      {
			streamlog_out(ERROR) << "Failed to extract a track: " << statusCodeException.ToString() << std::endl;
			streamlog_out( DEBUG3 ) << " failed track : " << *pTrack << std::endl ;
		      }
		    catch (EVENT::Exception &exception)
		      {
			streamlog_out(WARNING) << "Failed to extract a vertex: " << exception.what() << std::endl;
		      }
            }
            streamlog_out( DEBUG5 ) << "After treating collection : " << *iter<<" with "<<pTrackCollection->getNumberOfElements()<<" tracks, the track vector size is "<< m_trackVector.size()<< std::endl ;

        }
        catch (EVENT::Exception &exception)
        {
            streamlog_out(WARNING) << "Failed to extract track collection: " << *iter << ", " << exception.what() << std::endl;
        }
    }

    return pandora::STATUS_CODE_SUCCESS;
}

bool DDTrackCreatorALLEGRO::PassesQualityCuts(const EVENT::Track *const pTrack, const PandoraApi::Track::Parameters &trackParameters) const
{
    // First simple sanity checks
    if (trackParameters.m_trackStateAtCalorimeter.Get().GetPosition().GetMagnitude() < m_settings.m_minTrackECalDistanceFromIp){
        streamlog_out(WARNING) << " Dropping track! Distance at ECAL: " << trackParameters.m_trackStateAtCalorimeter.Get().GetPosition().GetMagnitude()<<std::endl;
	streamlog_out(DEBUG3)  << " track : " << *pTrack
			       << std::endl;
        return false;
    }


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
                               << " from " << pTrack->getTrackerHits().size()  << std::endl ;

	streamlog_out(DEBUG3)  << " track : " << *pTrack
			       << std::endl;
        return false;
    }

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DDTrackCreatorALLEGRO::DefineTrackPfoUsage(const EVENT::Track *const pTrack, PandoraApi::Track::Parameters &trackParameters) const
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

            const float zCutForNonVertexTracks(m_trackerInnerR * std::fabs(pZ / pT) + m_settings.m_zCutForNonVertexTracks);
            const bool passRzQualityCuts((zMin < zCutForNonVertexTracks) && (rInner < m_trackerInnerR + m_settings.m_maxBarrelTrackerInnerRDistance));

            const bool isV0(this->IsV0(pTrack));
            const bool isDaughter(this->IsDaughter(pTrack));

            // Decide whether track can be associated with a pandora cluster and used to form a charged PFO
            //if ((d0 < m_settings.m_d0TrackCut) && (z0 < m_settings.m_z0TrackCut) && (rInner < m_trackerInnerR + m_settings.m_maxBarrelTrackerInnerRDistance))
            // FIXME! AD: (rInner < m_trackerInnerR + m_settings.m_maxBarrelTrackerInnerRDistance) check is removed
            if ((d0 < m_settings.m_d0TrackCut) && (z0 < m_settings.m_z0TrackCut))
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
                    (rInner < m_trackerInnerR + m_settings.m_maxBarrelTrackerInnerRDistance))
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

void DDTrackCreatorALLEGRO::TrackReachesECAL(const EVENT::Track *const pTrack, PandoraApi::Track::Parameters &trackParameters) const
{
    // FIXME! AD: since currently we have do not have track hits -> check the radius of the trackState AtCalo.
    // currently check is done only for ECAL barrel but should check if track reaches the Endcap!
    pandora::CartesianVector posAtCalo(trackParameters.m_trackStateAtCalorimeter.Get().GetPosition());
    double radiusAtCalo = std::sqrt(posAtCalo.GetX()*posAtCalo.GetX() + posAtCalo.GetY()*posAtCalo.GetY());
    if(radiusAtCalo >= m_settings.m_eCalBarrelInnerR) trackParameters.m_reachesCalorimeter = true;
    return;


    // Calculate hit position information
    float hitZMin(std::numeric_limits<float>::max());
    float hitZMax(-std::numeric_limits<float>::max());
    float hitOuterR(-std::numeric_limits<float>::max());

    int nTrackerHits(0);
    int nEndcapDiskHits(0);
    int maxOccupiedEndcapDiskLayer(0);

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

        if ((r > m_trackerInnerR) && (r < m_trackerOuterR) && (std::fabs(z) <= m_trackerZmax))
        {
            nTrackerHits++;
            continue;
        }

        for (unsigned int j = 0; j < m_nEndcapDiskLayers; ++j)
        {
            if ((r > m_endcapDiskInnerRadii[j]) && (r < m_endcapDiskOuterRadii[j]) &&
                (std::fabs(z) - m_settings.m_reachesECalFtdZMaxDistance < m_endcapDiskZPositions[j]) &&
                (std::fabs(z) + m_settings.m_reachesECalFtdZMaxDistance > m_endcapDiskZPositions[j]))
            {
                if (static_cast<int>(j) > maxOccupiedEndcapDiskLayer)
                    maxOccupiedEndcapDiskLayer = static_cast<int>(j);

                nEndcapDiskHits++;
                break;
            }
        }
    }

    // Require sufficient hits in barrel or endcap trackers, then compare extremal hit positions with tracker dimensions
    if ((nTrackerHits >= m_settings.m_reachesECalNBarrelTrackerHits) || (nEndcapDiskHits >= m_settings.m_reachesECalNFtdHits))
    {
        if ((hitOuterR - m_trackerOuterR > m_settings.m_reachesECalBarrelTrackerOuterDistance) ||
            (std::fabs(hitZMax) - m_trackerZmax > m_settings.m_reachesECalBarrelTrackerZMaxDistance) ||
            (std::fabs(hitZMin) - m_trackerZmax > m_settings.m_reachesECalBarrelTrackerZMaxDistance) ||
            (maxOccupiedEndcapDiskLayer >= m_settings.m_reachesECalMinFtdLayer))
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

    if ((cosAngleAtDca > m_cosTracker) || (pT < m_settings.m_curvatureToMomentumFactor * m_settings.m_bField * m_trackerOuterR))
    {
        trackParameters.m_reachesCalorimeter = true;
        return;
    }

    trackParameters.m_reachesCalorimeter = false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DDTrackCreatorALLEGRO::GetTrackStatesAtCalo( EVENT::Track *track,
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
      // GM FIXME: for the moment, do not extrapolate also to endcap, and return
      return;
  } else { // if track state is in endcap we do not repeat track state calculation, because the barrel cannot be hit
      streamlog_out(DEBUG5) << "Original track state is at Endcap" << std::endl;
      pandora::InputTrackState pandoraTrackState;
      this->CopyTrackState( trackAtCalo, pandoraTrackState );
      trackParameters.m_trackStates.push_back( pandoraTrackState );
      return;
  }

  // GM FIXME: for the moment, do not extrapolate also to endcap, and return
  /*
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
  */
}

