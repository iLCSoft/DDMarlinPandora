/**
 *  @file   MarlinPandora/src/DDTrackCreatorILD.cc
 * 
 *  @brief  Implementation of the track creator class for ILD.
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

#include "DDTrackCreatorILD.h"
#include "Pandora/PdgTable.h"
#include "LCObjects/LCTrack.h"

#include <algorithm>
#include <cmath>
#include <limits>

#include "DD4hep/Detector.h"
#include "DD4hep/DD4hepUnits.h"
#include "DDRec/DetectorData.h"
#include "DD4hep/DetType.h"
#include "DD4hep/DetectorSelector.h"


DDTrackCreatorILD::DDTrackCreatorILD(const Settings &settings, const pandora::Pandora *const pPandora)
  :
  DDTrackCreatorBase(settings,pPandora),
  m_cosTpc( 0.f ),
  m_tpcInnerR( 0.f ),
  m_tpcOuterR( 0.f ),
  m_tpcMaxRow( 0 ),
  m_tpcZmax( 0.f ),
  m_tpcMembraneMaxZ( 0.f ),
  m_ftdInnerRadii( DoubleVector() ),
  m_ftdOuterRadii( DoubleVector() ),
  m_ftdZPositions( DoubleVector() ),
  m_nFtdLayers( 0 ),
  m_tanLambdaFtd( 0.f ),
  m_minEtdZPosition( 0.f ),
  m_minSetRadius( 0.f )

{
    
    
  m_nFtdLayers=0;
  m_ftdInnerRadii.clear();
  m_ftdOuterRadii.clear();
    
  //Instead of gear, loop over a provided list of forward (read: endcap) tracking detectors. For ILD this would be FTD
  ///FIXME: Should we use surfaces instead?
  dd4hep::Detector & mainDetector = dd4hep::Detector::getInstance();
  const std::vector< dd4hep::DetElement>& endcapDets = dd4hep::DetectorSelector(mainDetector).detectors(  ( dd4hep::DetType::TRACKER | dd4hep::DetType::ENDCAP )) ;
  for (std::vector< dd4hep::DetElement>::const_iterator iter = endcapDets.begin(), iterEnd = endcapDets.end();iter != iterEnd; ++iter){
    try
      {
	dd4hep::rec::ZDiskPetalsData * theExtension = 0;
  
        const dd4hep::DetElement& theDetector = *iter;
	theExtension = theDetector.extension<dd4hep::rec::ZDiskPetalsData>();
            
	unsigned int N = theExtension->layers.size();
            
	streamlog_out( DEBUG2 ) << " Filling FTD-like parameters from DD4hep for "<< theDetector.name() << "- n layers: " << N<< std::endl;

	for(unsigned int i = 0; i < N; ++i)
	  {
                
	    dd4hep::rec::ZDiskPetalsData::LayerLayout thisLayer  = theExtension->layers[i];

	    // Create a disk to represent even number petals front side
	    //FIXME! VERIFY THAT TIS MAKES SENSE!
	    m_ftdInnerRadii.push_back(thisLayer.distanceSensitive/dd4hep::mm);
	    m_ftdOuterRadii.push_back(thisLayer.distanceSensitive/dd4hep::mm+thisLayer.lengthSensitive/dd4hep::mm);

	    // Take the mean z position of the staggered petals
	    const double zpos(thisLayer.zPosition/dd4hep::mm);
	    m_ftdZPositions.push_back(zpos);
                
	    streamlog_out( DEBUG2 ) << "     layer " << i << " - mean z position = " << zpos << std::endl;
	  }

	m_nFtdLayers = m_ftdZPositions.size() ;
       
      } catch (std::runtime_error &exception){
            
      streamlog_out(WARNING) << "DDTrackCreatorILD exception during Forward Tracking Disk parameter construction for detector "<<const_cast<dd4hep::DetElement&>(*iter).name()<<" : " << exception.what() << std::endl;
    }
  }
    
  try
    {
      dd4hep::rec::FixedPadSizeTPCData * theExtension = 0;
      //Get the TPC, make sure not to get the vertex
      const std::vector< dd4hep::DetElement>& tpcDets= dd4hep::DetectorSelector(mainDetector).detectors(  ( dd4hep::DetType::TRACKER |  dd4hep::DetType::BARREL  | dd4hep::DetType::GASEOUS ), dd4hep::DetType::VERTEX) ;

      //There should only be one TPC
      theExtension = tpcDets[0].extension<dd4hep::rec::FixedPadSizeTPCData>();

      m_tpcInnerR = theExtension->rMin/dd4hep::mm ;
      m_tpcOuterR = theExtension->rMax/dd4hep::mm ;
      m_tpcMaxRow = theExtension->maxRow;
      m_tpcZmax   = theExtension->zHalf/dd4hep::mm ;
        
        
      //FIXME! Add to reco structrure and access
      m_tpcMembraneMaxZ = 10; 
      std::cout<<"WARNING! DO NOT MASK! HANDLE m_tpcMembraneMaxZ (currently hardcoded to 10)!"<<std::endl;
         
        
    } catch (std::runtime_error &exception){
    streamlog_out(WARNING) << "DDTrackCreatorILD exception during TPC parameter construction."<< std::endl;
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
  //FIXME: THINK OF A UNIVERSAL WAY TO HANDLE EXISTENCE OF ADDITIONAL DETECTORS

  streamlog_out(WARNING) << " ETDLayerZ or SETLayerRadius parameters Not being handled!" << std::endl
			 << "     -> both will be set to " << std::numeric_limits<float>::quiet_NaN() << std::endl;
  m_minEtdZPosition = std::numeric_limits<float>::quiet_NaN();
  m_minSetRadius = std::numeric_limits<float>::quiet_NaN();
    
  //     try
  //     {
  //         const DoubleVector &etdZPositions(marlin::Global::GEAR->getGearParameters("ETD").getDoubleVals("ETDLayerZ"));
  //         const DoubleVector &setInnerRadii(marlin::Global::GEAR->getGearParameters("SET").getDoubleVals("SETLayerRadius"));
  // 
  //         if (etdZPositions.empty() || setInnerRadii.empty())
  //             throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);
  // 
  //         m_minEtdZPosition = *(std::min_element(etdZPositions.begin(), etdZPositions.end()));
  //         m_minSetRadius = *(std::min_element(setInnerRadii.begin(), setInnerRadii.end()));
  //     }
  //     catch(gear::UnknownParameterException &)
  //     {
  //         streamlog_out(WARNING) << " ETDLayerZ or SETLayerRadius parameters missing from GEAR parameters!" << std::endl
  //                                << "     -> both will be set to " << std::numeric_limits<float>::quiet_NaN() << std::endl;
  // 
  //         //fg: Set them to NAN, so that they cannot be used to set   trackParameters.m_reachesCalorimeter = true;
  //         m_minEtdZPosition = std::numeric_limits<float>::quiet_NaN();
  //         m_minSetRadius = std::numeric_limits<float>::quiet_NaN();
  //     }
}

//------------------------------------------------------------------------------------------------------------------------------------------

DDTrackCreatorILD::~DDTrackCreatorILD()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode DDTrackCreatorILD::CreateTracks(EVENT::LCEvent *pLCEvent)
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
	      lc_content::LCTrackParameters trackParameters;
	      trackParameters.m_d0 = pTrack->getD0();
	      trackParameters.m_z0 = pTrack->getZ0();
	      trackParameters.m_pParentAddress = pTrack;

	      // By default, assume tracks are charged pions
	      const float signedCurvature(pTrack->getOmega());
	      trackParameters.m_particleId = (signedCurvature > 0) ? pandora::PI_PLUS : pandora::PI_MINUS;
	      trackParameters.m_mass = pandora::PdgTable::GetParticleMass(pandora::PI_PLUS);

	      // Use particle id information from V0 and Kink finders
	      TrackToPidMap::const_iterator trIter = m_trackToPidMap.find(pTrack);

	      if(trIter != m_trackToPidMap.end()) {
		trackParameters.m_particleId = trIter->second;
		trackParameters.m_mass = pandora::PdgTable::GetParticleMass(trIter->second);
	      }

	      if (std::numeric_limits<float>::epsilon() < std::fabs(signedCurvature))
		trackParameters.m_charge = static_cast<int>(signedCurvature / std::fabs(signedCurvature));

	      try //fg: include the next calls in the try block to catch tracks that are yet not fitted properly as ERROR and not exceptions
		{
		  this->GetTrackStates(pTrack, trackParameters);
		  this->TrackReachesECAL(pTrack, trackParameters);
		  this->GetTrackStatesAtCalo(pTrack, trackParameters);
		  this->DefineTrackPfoUsage(pTrack, trackParameters);

		  PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::Track::Create(m_pandora, trackParameters, *m_lcTrackFactory));
		  m_trackVector.push_back(pTrack);
		}
	      catch (pandora::StatusCodeException &statusCodeException)
		{
		  streamlog_out(ERROR) << "Failed to extract a track: " << statusCodeException.ToString() << std::endl;
			
		  streamlog_out( DEBUG5 ) << " failed track : " << *pTrack << std::endl ;
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

bool DDTrackCreatorILD::PassesQualityCuts(const EVENT::Track *const pTrack, const PandoraApi::Track::Parameters &trackParameters) const
{
  // First simple sanity checks
  if (trackParameters.m_trackStateAtCalorimeter.Get().GetPosition().GetMagnitude() < m_settings.m_minTrackECalDistanceFromIp)
    return false;

  if (std::fabs(pTrack->getOmega()) < std::numeric_limits<float>::epsilon())
    {
      streamlog_out(ERROR) << "Track has Omega = 0 " << std::endl;
      return false;
    }


  if( pTrack->getNdf() < 0 ){
    streamlog_out(ERROR) << "Track is unconstrained - ndf = " <<  pTrack->getNdf()  << std::endl;
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

      streamlog_out(DEBUG5)  << " track : " << *pTrack
			     << std::endl;
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

      if ((std::numeric_limits<float>::epsilon() > std::fabs(pT)) || (std::numeric_limits<float>::epsilon() > std::fabs(pZ)) || (rInnermostHit == m_tpcOuterR))
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
      if (std::fabs(pZ) / momentumAtDca.GetMagnitude() < m_tpcMembraneMaxZ / m_tpcInnerR)
	nExpectedTpcHits = 0;

      const int nTpcHits(this->GetNTpcHits(pTrack));
      const int nFtdHits(this->GetNFtdHits(pTrack));

      const int minTpcHits = static_cast<int>(nExpectedTpcHits * m_settings.m_minBarrelTrackerHitFractionOfExpected);

      if ((nTpcHits < minTpcHits) && (nFtdHits < m_settings.m_minFtdHitsForBarrelTrackerHitFraction))
        {
	  streamlog_out(WARNING) << " Dropping track : " << momentumAtDca.GetMagnitude() << " Number of TPC hits = " << nTpcHits
				 << " < " << minTpcHits << " nftd = " << nFtdHits << std::endl;

	  streamlog_out(DEBUG5)  << " track : " << *pTrack
				 << std::endl;
	  return false;
        }
    }

  return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DDTrackCreatorILD::DefineTrackPfoUsage(const EVENT::Track *const pTrack, PandoraApi::Track::Parameters &trackParameters) const
{
  bool canFormPfo(false);
  bool canFormClusterlessPfo(false);

  if (trackParameters.m_reachesCalorimeter.Get() && !this->IsParent(pTrack))
    {
      const float d0(std::fabs(pTrack->getD0())), z0(std::fabs(pTrack->getZ0()));

      streamlog_out(DEBUG3)  << " -- DefineTrackPfoUsage called for track : " << UTIL::lcshort( pTrack )
			     << std::endl;

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
	  const bool passRzQualityCuts((zMin < zCutForNonVertexTracks) && (rInner < m_tpcInnerR + m_settings.m_maxBarrelTrackerInnerRDistance));

	  const bool isV0(this->IsV0(pTrack));
	  const bool isDaughter(this->IsDaughter(pTrack));

	  streamlog_out(DEBUG3)  << " -- track passed quality cuts and has : " 
				 << " passRzQualityCuts " << passRzQualityCuts
				 << " isV0 " << isV0
				 << " isDaughter " << isDaughter << std::endl ;

	  // Decide whether track can be associated with a pandora cluster and used to form a charged PFO
	  if ((d0 < m_settings.m_d0TrackCut) && (z0 < m_settings.m_z0TrackCut) && (rInner < m_tpcInnerR + m_settings.m_maxBarrelTrackerInnerRDistance))
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
		  (rInner < m_tpcInnerR + m_settings.m_maxBarrelTrackerInnerRDistance))
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
      
      streamlog_out(DEBUG3)  << " -- track canFormPfo = " << canFormPfo << " -  canFormClusterlessPfo = " << canFormClusterlessPfo << std::endl;

    }

  trackParameters.m_canFormPfo = canFormPfo;
  trackParameters.m_canFormClusterlessPfo = canFormClusterlessPfo;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DDTrackCreatorILD::TrackReachesECAL(const EVENT::Track *const pTrack, PandoraApi::Track::Parameters &trackParameters) const
{
  //fg: return true  for now - there are quality checks in DefineTrackPfoUsage() ...
  trackParameters.m_reachesCalorimeter = true;
  return;
  
  // at a later stage we could simply check if there is a valid track state at the calorimeter:
  // ...
  
  
  
  
  // Calculate hit position information
  float hitZMin(std::numeric_limits<float>::max());
  float hitZMax(-std::numeric_limits<float>::max());
  float hitOuterR(-std::numeric_limits<float>::max());
  
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

	      break;
            }
        }
    }

    const int nTpcHits(this->GetNTpcHits(pTrack));
    const int nFtdHits(this->GetNFtdHits(pTrack));
    
  // Look to see if there are hits in etd or set, implying track has reached edge of ecal
  if ((hitOuterR > m_minSetRadius) || (hitZMax > m_minEtdZPosition))
    {
      trackParameters.m_reachesCalorimeter = true;
      return;
    }

  // Require sufficient hits in tpc or ftd, then compare extremal hit positions with tracker dimensions
  if ((nTpcHits >= m_settings.m_reachesECalNBarrelTrackerHits) || (nFtdHits >= m_settings.m_reachesECalNFtdHits))
    {
      if ((hitOuterR - m_tpcOuterR > m_settings.m_reachesECalBarrelTrackerOuterDistance) ||
	  (std::fabs(hitZMax) - m_tpcZmax > m_settings.m_reachesECalBarrelTrackerZMaxDistance) ||
	  (std::fabs(hitZMin) - m_tpcZmax > m_settings.m_reachesECalBarrelTrackerZMaxDistance) ||
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

  if ((cosAngleAtDca > m_cosTpc) || (pT < m_settings.m_curvatureToMomentumFactor * m_settings.m_bField * m_tpcOuterR))
    {
      trackParameters.m_reachesCalorimeter = true;
      return;
    }

  trackParameters.m_reachesCalorimeter = false;
}

int DDTrackCreatorILD::GetNTpcHits(const EVENT::Track *const pTrack) const
{
   // ATTN
   //According to FG: [ 2 * lcio::ILDDetID::TPC - 2 ] is the first number and it is supposed to
   //be the number of hits in the fit and this is what should be used !
   // at least for DD4hep/DDSim

    // ---- use hitsInFit :
    return pTrack->getSubdetectorHitNumbers()[ 2 * lcio::ILDDetID::TPC - 2 ];
}

//------------------------------------------------------------------------------------------------------------------------------------------

int DDTrackCreatorILD::GetNFtdHits(const EVENT::Track *const pTrack) const
{
    // ATTN
    //see above
    // ---- use hitsInFit :
    return pTrack->getSubdetectorHitNumbers()[ 2 * lcio::ILDDetID::FTD - 2 ];
}

//
