/**
 *  @file   MarlinPandora/src/DDTrackCreatorCLIC.cc
 * 
 *  @brief  Implementation of the track creator class for a CLIC all silicon tracker.
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
#include "DDTrackCreatorCLIC.h"
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

DDTrackCreatorCLIC::DDTrackCreatorCLIC(const Settings &settings, const pandora::Pandora *const pPandora) : DDTrackCreatorBase(settings,pPandora)
{
    
    m_trackerInnerR = getTrackingRegionExtent()[0];
    m_trackerOuterR = getTrackingRegionExtent()[1];
    m_trackerZmax = getTrackingRegionExtent()[2];

    
    ///FIXME: Probably need to be something relating to last disk inner radius
    m_cosTracker = m_trackerZmax / std::sqrt(m_trackerZmax * m_trackerZmax + m_trackerInnerR * m_trackerInnerR);

    
    
    m_barrelTrackerLayers = 0;
    for (std::vector<std::string>::const_iterator iter = m_settings.m_barrelTrackerNames.begin(), iterEnd = m_settings.m_barrelTrackerNames.end();iter != iterEnd; ++iter){
        try
        {
            DD4hep::DDRec::ZPlanarData * theExtension = 0;
  
            DD4hep::Geometry::LCDD & lcdd = DD4hep::Geometry::LCDD::getInstance();
            DD4hep::Geometry::DetElement theDetector = lcdd.detector(*iter);
            theExtension = theDetector.extension<DD4hep::DDRec::ZPlanarData>();
            
            unsigned int N = theExtension->layers.size();
            m_barrelTrackerLayers=m_barrelTrackerLayers+N;
            
            streamlog_out( DEBUG2 ) << " Adding layers for barrel tracker from DD4hep for "<< *iter<< "- n layers: " << N<< " sum up to now: "<<m_barrelTrackerLayers<<std::endl;
        } catch (std::runtime_error &exception){
            
            streamlog_out(WARNING) << "DDTrackCreatorCLIC exception during Barrel Tracker layer sum for "<<*iter<<" : " << exception.what() << std::endl;
        }
    }
    
    m_nFtdLayers=0;
    m_ftdInnerRadii.clear();
    m_ftdOuterRadii.clear();
    
    //Instead of gear, loop over a provided list of forward (read: endcap) tracking detectors. For ILD this would be FTD
    ///FIXME: Should we use surfaces instead?
    for (std::vector<std::string>::const_iterator iter = m_settings.m_endcapTrackerNames.begin(), iterEnd = m_settings.m_endcapTrackerNames.end();iter != iterEnd; ++iter){
        try
        {
            DD4hep::DDRec::ZDiskPetalsData * theExtension = 0;
  
            DD4hep::Geometry::LCDD & lcdd = DD4hep::Geometry::LCDD::getInstance();
            DD4hep::Geometry::DetElement theDetector = lcdd.detector(*iter);
            theExtension = theDetector.extension<DD4hep::DDRec::ZDiskPetalsData>();
            
            unsigned int N = theExtension->layers.size();
            
            streamlog_out( DEBUG2 ) << " Filling FTD-like parameters from DD4hep for "<< *iter<< "- n layers: " << N<< std::endl;

            for(unsigned int i = 0; i < N; ++i)
            {
                
                DD4hep::DDRec::ZDiskPetalsData::LayerLayout thisLayer  = theExtension->layers[i];

                // Create a disk to represent even number petals front side
                //FIXME! VERIFY THAT TIS MAKES SENSE!
                m_ftdInnerRadii.push_back(thisLayer.distanceSensitive);
                m_ftdOuterRadii.push_back(thisLayer.distanceSensitive+thisLayer.lengthSensitive);

                // Take the mean z position of the staggered petals
                const double zpos(thisLayer.zPosition);
                m_ftdZPositions.push_back(zpos);
                
                streamlog_out( DEBUG2 ) << "     layer " << i << " - mean z position = " << zpos << std::endl;
            }

            m_nFtdLayers = m_ftdZPositions.size() ;
       
        } catch (std::runtime_error &exception){
            
            streamlog_out(WARNING) << "DDTrackCreatorCLIC exception during Forward Tracking Disk parameter construction for detector "<<*iter<<" : " << exception.what() << std::endl;
        }
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

		    try
		      {
			PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::Track::Create(*m_pPandora, trackParameters));
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

bool DDTrackCreatorCLIC::PassesQualityCuts(const EVENT::Track *const pTrack, const PandoraApi::Track::Parameters &trackParameters) const
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
                               << " from " << pTrack->getTrackerHits().size()  << std::endl ;

	streamlog_out(DEBUG5)  << " track : " << *pTrack
			       << std::endl;
        return false;
    }

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
        //FIXME! NEED TO COME UP WITH CONSISTENT CONVENTION FOR INNER/OUTER TRACKER
        const int nTrackerHits = hitsBySubdetector[ 2 * lcio::ILDDetID::TPC - 1 ];
        const int nFtdHits = hitsBySubdetector[ 2 * lcio::ILDDetID::FTD - 1 ];

        const int minTrackerHits = static_cast<int>(nExpectedTrackerHits * m_settings.m_minBarrelTrackerHitFractionOfExpected);

        if ((nTrackerHits < minTrackerHits) && (nFtdHits < m_settings.m_minFtdHitsForBarrelTrackerHitFraction))
        {
            streamlog_out(WARNING) << " Dropping track : " << momentumAtDca.GetMagnitude() << " Number of Tracker hits = " << nTrackerHits
                                   << " < " << minTrackerHits << " nftd = " << nFtdHits << std::endl;

	    streamlog_out(DEBUG5)  << " track : " << *pTrack
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

        if ((r > m_trackerInnerR) && (r < m_trackerOuterR) && (std::fabs(z) <= m_trackerZmax))
        {
            nTrackerHits++;
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

    // Require sufficient hits in tpc or ftd, then compare extremal hit positions with tracker dimensions
    if ((nTrackerHits >= m_settings.m_reachesECalNBarrelTrackerHits) || (nFtdHits >= m_settings.m_reachesECalNFtdHits))
    {
        if ((hitOuterR - m_trackerOuterR > m_settings.m_reachesECalBarrelTrackerOuterDistance) ||
            (std::fabs(hitZMax) - m_trackerZmax > m_settings.m_reachesECalBarrelTrackerZMaxDistance) ||
            (std::fabs(hitZMin) - m_trackerZmax > m_settings.m_reachesECalBarrelTrackerZMaxDistance) ||
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

    if ((cosAngleAtDca > m_cosTracker) || (pT < m_settings.m_curvatureToMomentumFactor * m_settings.m_bField * m_trackerOuterR))
    {
        trackParameters.m_reachesCalorimeter = true;
        return;
    }

    trackParameters.m_reachesCalorimeter = false;
}