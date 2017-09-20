/**
 *  @file   DDMarlinPandora/src/DDTrackCreatorCLIC.cc
 * 
 *  @brief  Implementation of the track creator class for a CLIC all silicon tracker.
 * 
 *  $Log: $
 */

#include "DDTrackCreatorCLIC.h"

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

DDTrackCreatorCLIC::DDTrackCreatorCLIC(const Settings &settings, const pandora::Pandora *const pPandora)
  : DDTrackCreatorBase(settings,pPandora),
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
            
            streamlog_out(WARNING) << "DDTrackCreatorCLIC exception during Barrel Tracker layer sum for "<<const_cast<dd4hep::DetElement&>(*iter).name()<<" : " << exception.what() << std::endl;
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
            
            streamlog_out(WARNING) << "DDTrackCreatorCLIC exception during Forward Tracking Disk parameter construction for detector "<<const_cast<dd4hep::DetElement&>(*iter).name()<<" : " << exception.what() << std::endl;
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

}

//------------------------------------------------------------------------------------------------------------------------------------------

DDTrackCreatorCLIC::~DDTrackCreatorCLIC()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode DDTrackCreatorCLIC::CreateTracks(EVENT::LCEvent *pLCEvent)
{

    for (StringVector::const_iterator iter = m_settings.m_trackCollections.begin(), iterEnd = m_settings.m_trackCollections.end();
        iter != iterEnd; ++iter)
    {
        try
        {
            const EVENT::LCCollection *pTrackCollection = pLCEvent->getCollection(*iter);
            
            
            ///FIXME: Should really move to using surfaces
            for (int i = 0, iMax = pTrackCollection->getNumberOfElements(); i < iMax; ++i)
            {
                    EVENT::Track *pTrack = dynamic_cast<Track*>(pTrackCollection->getElementAt(i));

                    if (NULL == pTrack)
                        throw EVENT::Exception("Collection type mismatch");

		    streamlog_out(DEBUG0)<<" Warning! Ignoring expected number of hits and other hit number cuts. Should eventually change!"<<std::endl;

//                     int minTrackHits = m_settings.m_minTrackHits;
//                     const float tanLambda(std::fabs(pTrack->getTanLambda()));
// 
//                     if (tanLambda > m_tanLambdaEndcapDisk)
//                     {
//                         int expectedEndcapDiskHits(0);
// 
//                         for (unsigned int iEndcapDiskLayer = 0; iEndcapDiskLayer < m_nEndcapDiskLayers; ++iEndcapDiskLayer)
//                         {
//                             
//                             //FIXME: Does not take into account spiral endcap
//                             if ((tanLambda > m_endcapDiskZPositions[iEndcapDiskLayer] / m_endcapDiskOuterRadii[iEndcapDiskLayer]) &&
//                                 (tanLambda < m_endcapDiskZPositions[iEndcapDiskLayer] / m_endcapDiskInnerRadii[iEndcapDiskLayer]))
//                             {
//                                 expectedEndcapDiskHits++;
//                             }
//                         }
// 
//                         minTrackHits = std::max(m_settings.m_minEndcapDiskTrackHits, expectedEndcapDiskHits);
//                         
//                         streamlog_out(DEBUG0)<<"XXX minTrackHits: "<<minTrackHits<<" m_minEndcapDiskTrackHits: "<<m_settings.m_minEndcapDiskTrackHits<<" expectedEndcapDiskHits: "<<expectedEndcapDiskHits<<std::endl;
//                     }
// 
//                     const int nTrackHits(static_cast<int>(pTrack->getTrackerHits().size()));
// 
//                     if ((nTrackHits < minTrackHits) || (nTrackHits > m_settings.m_maxTrackHits)){
//                         streamlog_out(DEBUG0)<<"XXX Dropping track due to nTrackHits= "<<nTrackHits<<" . minTrackHits: "<<minTrackHits<<" maxTrackHits: "<<m_settings.m_maxTrackHits<<std::endl;
//                         continue;
//                     }

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
			//FIXME: Should consider adding a try-catch block to check
			//against cases where a track has invalid parameters like
			//very small omega. Especially for omega<epsilon=1e-7,
			//pandora throws a pandora::StatusCodeException
			//&statusCodeException which is not caught.
			//FIXED: if an exception happens we ignore the track

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

bool DDTrackCreatorCLIC::PassesQualityCuts(const EVENT::Track *const pTrack, const PandoraApi::Track::Parameters &trackParameters) const
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

    streamlog_out( DEBUG3 ) << " TEMPORARILY ACCEPT TRACK WITHOUT CUTS (should change!)" << *pTrack << std::endl ;
    return true;
    
    // Require reasonable number of Tracker hits 
    if (momentumAtDca.GetMagnitude() > m_settings.m_minMomentumForTrackHitChecks)
    {
        const float pX(fabs(momentumAtDca.GetX()));
        const float pY(fabs(momentumAtDca.GetY()));
        const float pZ(fabs(momentumAtDca.GetZ()));
        const float pT(std::sqrt(pX * pX + pY * pY));
        const float rInnermostHit(pTrack->getRadiusOfInnermostHit());

        if ((0.f == pT) || (0.f == pZ) || (rInnermostHit == m_trackerOuterR))
        {
            streamlog_out(ERROR) << "Invalid track parameter, pT " << pT << ", pZ " << pZ << ", rInnermostHit " << rInnermostHit << std::endl;
            return false;
        }

        float nExpectedTrackerHits(0.);

        if (pZ < m_trackerZmax / m_trackerOuterR * pT)
        {
            const float innerExpectedHitRadius(std::max(m_trackerInnerR, rInnermostHit));
            const float frac((m_trackerOuterR - innerExpectedHitRadius) / (m_trackerOuterR - m_trackerInnerR));
            nExpectedTrackerHits = m_barrelTrackerLayers * frac;
        }

        if ((pZ <= m_trackerZmax / m_trackerInnerR * pT) && (pZ >= m_trackerZmax / m_trackerOuterR * pT))
        {
            const float innerExpectedHitRadius(std::max(m_trackerInnerR, rInnermostHit));
            const float frac((m_trackerZmax * pT / pZ - innerExpectedHitRadius) / (m_trackerOuterR - innerExpectedHitRadius));
            nExpectedTrackerHits = frac * m_barrelTrackerLayers;
        }



        const EVENT::IntVec &hitsBySubdetector(pTrack->getSubdetectorHitNumbers());
 
        // ---- use hitsInFit :
        //Use flags to access detectors and then access their ids to get hits in fit by ID
        //(new way of storing hits in MarlinTrkUtils)
        
        //Initialize hits to 0
        int nBarrelTrackerHits = 0;
        dd4hep::Detector & mainDetector = dd4hep::Detector::getInstance();
        const std::vector< dd4hep::DetElement>& barrelDets = dd4hep::DetectorSelector(mainDetector).detectors(  ( dd4hep::DetType::TRACKER | dd4hep::DetType::BARREL )) ;
        for (std::vector< dd4hep::DetElement>::const_iterator iter = barrelDets.begin(), iterEnd = barrelDets.end();iter != iterEnd; ++iter){
            const dd4hep::DetElement & theDetector = *iter;
            int detId = theDetector.id();
            nBarrelTrackerHits+=hitsBySubdetector[2*detId-2]; //Offset is 2 for hits in fit
        }
        
        int nEndcapTrackerHits = 0;
        const std::vector< dd4hep::DetElement>& endcapDets = dd4hep::DetectorSelector(mainDetector).detectors(  ( dd4hep::DetType::TRACKER | dd4hep::DetType::ENDCAP )) ;
        for (std::vector< dd4hep::DetElement>::const_iterator iter = endcapDets.begin(), iterEnd = endcapDets.end();iter != iterEnd; ++iter){
            const dd4hep::DetElement& theDetector = *iter;
            int detId = theDetector.id();
            nEndcapTrackerHits +=hitsBySubdetector[2*detId-2]; //Offset is 2 for hits in fit
        }

        const int minTrackerHits = static_cast<int>(nExpectedTrackerHits * m_settings.m_minBarrelTrackerHitFractionOfExpected);

        if ((nBarrelTrackerHits < minTrackerHits) && (nEndcapTrackerHits < m_settings.m_minFtdHitsForBarrelTrackerHitFraction))
        {
            streamlog_out(WARNING) << " Dropping track : " << momentumAtDca.GetMagnitude() << " Number of Tracker hits = " << nBarrelTrackerHits
                                   << " < " << minTrackerHits << " nendcapDisk = " << nEndcapTrackerHits << std::endl;

	    streamlog_out(DEBUG3)  << " track : " << *pTrack
				   << std::endl;
            return false;
        }
    }

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DDTrackCreatorCLIC::DefineTrackPfoUsage(const EVENT::Track *const pTrack, PandoraApi::Track::Parameters &trackParameters) const
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
            if ((d0 < m_settings.m_d0TrackCut) && (z0 < m_settings.m_z0TrackCut) && (rInner < m_trackerInnerR + m_settings.m_maxBarrelTrackerInnerRDistance))
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

void DDTrackCreatorCLIC::TrackReachesECAL(const EVENT::Track *const pTrack, PandoraApi::Track::Parameters &trackParameters) const
{
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
