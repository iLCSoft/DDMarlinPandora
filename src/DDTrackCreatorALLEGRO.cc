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
#include "UTIL/Operators.h"

#include "Pandora/PdgTable.h"
#include "LCObjects/LCTrack.h"

#include "DD4hep/Detector.h"
#include "DD4hep/DD4hepUnits.h"
#include "DDRec/DetectorData.h"
#include "DDRec/DCH_info.h"
#include "DD4hep/DetType.h"
#include "DD4hep/DetectorSelector.h"

#include <algorithm>
#include <cmath>
#include <limits>

//forward declarations. See in DDPandoraPFANewProcessor.cc
std::vector<double> getTrackingRegionExtent();

DDTrackCreatorALLEGRO::DDTrackCreatorALLEGRO(const Settings &settings, const pandora::Pandora *const pPandora)
  : DDTrackCreatorBase(settings, pPandora),
    m_dchInnerR( 0.f ),
    m_dchOuterR( 0.f ),
    m_dchOuterZ( 0.f ),
    m_cosDch( 0.f ),
    m_dchNLayers( 0 ),
    m_wrapperBarrelInnerR( 0.f ),
    m_wrapperBarrelOuterR( 0.f ),
    m_wrapperBarrelOuterZ( 0.f ),
    m_wrapperEndcapInnerR( 0.f ),
    m_wrapperEndcapOuterR( 0.f ),
    m_wrapperEndcapInnerZ( 0.f ),
    m_wrapperEndcapOuterZ( 0.f ),
    m_wrapperBarrelNLayers( 0 ),
    m_wrapperEndcapNLayers( 0 )
    /*
    m_barrelWrapperRPositions( DoubleVector() ),
    m_barrelWrapperOuterZ( DoubleVector() ),
    m_endcapWrapperInnerR( DoubleVector() ),
    m_endcapWrapperOuterR( DoubleVector() ),
    m_endcapWrapperZPositions( DoubleVector() ),
    m_nWrapperBarrelLayers( 0 ),
    m_nWrapperEndcapLayers( 0 ),
    m_nTrackerLayers( 0 ),
    m_tanLambdaEndcapDisk( 0.f )
*/
{

    dd4hep::Detector & mainDetector = dd4hep::Detector::getInstance();

    // Get DCH parameters
    try {
        const std::vector< dd4hep::DetElement>& dchDets= dd4hep::DetectorSelector(mainDetector).detectors(  ( dd4hep::DetType::TRACKER |  dd4hep::DetType::BARREL  | dd4hep::DetType::GASEOUS ), dd4hep::DetType::VERTEX) ;
        auto dchExtension = dchDets[0].extension<dd4hep::rec::DCH_info>();
        m_dchInnerR = dchExtension->rin/dd4hep::mm ;
        m_dchOuterR = dchExtension->rout/dd4hep::mm ;
        m_dchOuterZ = dchExtension->Lhalf/dd4hep::mm ;
        m_dchNLayers = dchExtension->nlayers;
    }
    catch (...) {
        std::cout << "Failed to retrieve DCH information" << std::endl;
    }
    if ((std::fabs(m_dchOuterZ) < std::numeric_limits<float>::epsilon()) ||
        (std::fabs(m_dchInnerR) < std::numeric_limits<float>::epsilon()) ||
        (std::fabs(m_dchOuterR - m_dchInnerR) < std::numeric_limits<float>::epsilon()))
    {
        throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);
    }
    m_cosDch = m_dchOuterZ / std::sqrt(m_dchOuterZ * m_dchOuterZ + m_dchInnerR * m_dchInnerR);



    // Get wrapper parameters
    // FIXME!

    // const std::vector< dd4hep::DetElement>& wrapperBarrelDets = dd4hep::DetectorSelector(mainDetector).detectors(  ( dd4hep::DetType::TRACKER | dd4hep::DetType::PIXEL | dd4hep::DetType::BARREL ), dd4hep::DetType::VERTEX) ;
    // d4hep::rec::ZPlanarData * wrBExtension = wrapperBarrelDets[0].extension<dd4hep::rec::ZPlanarData>();
    // m_wrBNLayers = wrBExtension->layers.size();
    // ..
    // const std::vector< dd4hep::DetElement>& wrapperEndcapDets = dd4hep::DetectorSelector(mainDetector).detectors(  ( dd4hep::DetType::TRACKER | dd4hep::DetType::PIXEL | dd4hep::DetType::ENDCAP ), dd4hep::DetType::VERTEX) ;
    // d4hep::rec::ZDiskPetalsData * wrEExtension = wrapperEndcapDets[0].extension<dd4hep::rec::ZDiskPetalsData>();
    // m_wrENLayers = wrEExtension->layers.size();
    // ..

    // GM: there is no reconstruction data saved for the wrapper yet...
    // so I just define its envelope based on the limits of the drift chamber and of the full tracking volume
    // full tracking volume
    float trackerInnerR = getTrackingRegionExtent()[0];
    float trackerOuterR = getTrackingRegionExtent()[1];
    float trackerOuterZ = getTrackingRegionExtent()[2];
    m_wrapperBarrelInnerR = m_dchOuterR + 1.0;
    m_wrapperBarrelOuterR = trackerOuterR;
    m_wrapperBarrelOuterZ = m_dchOuterZ;
    m_wrapperEndcapInnerR = m_dchInnerR;
    m_wrapperEndcapOuterR = trackerOuterR;
    m_wrapperEndcapInnerZ = m_dchOuterZ + 1.0;
    m_wrapperEndcapOuterZ = trackerOuterZ;
    m_wrapperBarrelNLayers = 2;
    m_wrapperEndcapNLayers = 2;
    m_cosWrapper = m_wrapperEndcapOuterZ / std::sqrt(m_wrapperEndcapOuterZ * m_wrapperEndcapOuterZ + m_wrapperEndcapInnerR * m_wrapperEndcapInnerR);

    // streamlog_out(DEBUG0)
    streamlog_out(MESSAGE)
        << "DDTrackCreatorALLEGRO DEBUG: " << std::endl
                          << "DCH rIn, rOut, zOut (mm): " << m_dchInnerR << " , " << m_dchOuterR << " , " << m_dchOuterZ << std::endl
                          << "Wrapper barrel rIn, rOut, zOut (mm): " << m_wrapperBarrelInnerR << " , " << m_wrapperBarrelOuterR << " , " << m_wrapperBarrelOuterZ << std::endl
                          << "Wrapper endcap rIn, rOut, zIn, zOut (mm): " << m_wrapperEndcapInnerR << " , " << m_wrapperEndcapOuterR << " , " << m_wrapperEndcapInnerZ << " , " << m_wrapperEndcapOuterZ << std::endl
                          << "number of layers in DCH, wrapper barrel, wrapper endcap: " << m_dchNLayers << " , " << m_wrapperBarrelNLayers << " , " << m_wrapperEndcapNLayers << std::endl;
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

                // GM note that anyway we have cuts on number of hits in the PassQualityCuts method.
                // maybe here just a very loose cut on the number of hits?
                int minTrackHits = m_settings.m_minTrackHits;
                int maxTrackHits = m_settings.m_maxTrackHits;
                const int nTrackHits(this->GetNVertexHits(pTrack) + this->GetNDchHits(pTrack) + this->GetNSiWrapperHits(pTrack));
                if ((nTrackHits < minTrackHits) || (nTrackHits > maxTrackHits)) {
                    streamlog_out(WARNING) << " Dropping track : " << " Number of hits = " << nTrackHits
                                           << " is not in [ " << minTrackHits << " , " << maxTrackHits << " ] range" << std::endl;
                    continue;
                }

                // Proceed to create the pandora track
                lc_content::LCTrackParameters trackParameters;
                trackParameters.m_d0 = pTrack->getD0();
                trackParameters.m_z0 = pTrack->getZ0();
                trackParameters.m_pParentAddress = pTrack;

                // By default, assume tracks are charged pions
                const float signedCurvature(pTrack->getOmega());
                trackParameters.m_particleId = (signedCurvature > 0) ? pandora::PI_PLUS : pandora::PI_MINUS;
                trackParameters.m_mass = pandora::PdgTable::GetParticleMass(pandora::PI_PLUS);

                // Use particle id information from V0 and Kink finders (if any)
                TrackToPidMap::const_iterator trackPIDiter = m_trackToPidMap.find(pTrack);
                if(trackPIDiter != m_trackToPidMap.end())
                {
                    trackParameters.m_particleId = trackPIDiter->second;
                    trackParameters.m_mass = pandora::PdgTable::GetParticleMass(trackPIDiter->second);
                }

                // Set charge if curvature different from zero
                if (0.f != signedCurvature)
                    trackParameters.m_charge = static_cast<int>(signedCurvature / std::fabs(signedCurvature));

                try { // include the next calls in the try block to catch tracks that are yet not fitted properly as ERROR and not exceptions
                    // retrieve track position at DCA, first and last hits, and at calorimeter
                    this->GetTrackStates(pTrack, trackParameters);
                    // determine if track reaches calorimeter
                    this->TrackReachesECAL(pTrack, trackParameters);
                    // calculate possible additional states at calorimeter
                    this->GetTrackStatesAtCalo(pTrack, trackParameters);
                    // decide if tracks can form a PFO with/without matched cluster
                    this->DefineTrackPfoUsage(pTrack, trackParameters);
                    // create Pandora track and add it to list of tracks
                    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::Track::Create(m_pandora, trackParameters, *m_lcTrackFactory));
                    m_trackVector.push_back(pTrack);
                }
                catch (pandora::StatusCodeException &statusCodeException) {
                    streamlog_out(ERROR) << "Failed to extract a track: " << statusCodeException.ToString() << std::endl;
                    std::cout << " failed track : " << *pTrack << std::endl ;
                    streamlog_out( DEBUG3 ) << " failed track : " << *pTrack << std::endl ;
                }
                catch (EVENT::Exception &exception) {
                    streamlog_out(WARNING) << "Failed to extract a vertex: " << exception.what() << std::endl;
                }
            }
            //streamlog_out( DEBUG5 )
            streamlog_out( MESSAGE )
                << "After treating collection : " << *iter<<" with "<<pTrackCollection->getNumberOfElements()<<" tracks, the track vector size is "<< m_trackVector.size()<< std::endl ;
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
    // - distance at ECAL from IP greater than minimum threshold
    if (trackParameters.m_trackStateAtCalorimeter.Get().GetPosition().GetMagnitude() < m_settings.m_minTrackECalDistanceFromIp){
        streamlog_out(WARNING) << " Dropping track! Distance at ECAL: " << trackParameters.m_trackStateAtCalorimeter.Get().GetPosition().GetMagnitude()<<std::endl;
        streamlog_out(DEBUG3)  << " track : " << *pTrack
                               << std::endl;
        return false;
    }
    // - non-zero curvature
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

    // Require reasonable number of DCH hits
    if (momentumAtDca.GetMagnitude() > m_settings.m_minMomentumForTrackHitChecks)
    {
        const float pX(fabs(momentumAtDca.GetX()));
        const float pY(fabs(momentumAtDca.GetY()));
        const float pZ(fabs(momentumAtDca.GetZ()));
        const float pT(std::sqrt(pX * pX + pY * pY));
        const float rInnermostHit(pTrack->getRadiusOfInnermostHit());

        // reject track with zero pT or pZ (why for pZ?) or innermost hit beyond DCH
        if ((std::numeric_limits<float>::epsilon() > std::fabs(pT)) || (std::numeric_limits<float>::epsilon() > std::fabs(pZ)) || (rInnermostHit >= m_dchOuterR))
        {
            streamlog_out(ERROR) << "Invalid track parameter, pT " << pT << ", pZ " << pZ << ", rInnermostHit " << rInnermostHit << std::endl;
            return false;
        }

        // calculate number of expected DCH hits based on track direction (pT, pZ), DCH sides and position of innermost track hit
        float nExpectedDchHits(0.);
        // case where a projective track from IP reaches the outer side of the DCH cylinder
        if (pZ < m_dchOuterZ / m_dchOuterR * pT)
        {
            const float innerExpectedHitRadius(std::max(m_dchInnerR, rInnermostHit));
            const float frac((m_dchOuterR - innerExpectedHitRadius) / (m_dchOuterR - m_dchInnerR));
            nExpectedDchHits = m_dchNLayers * frac;
        }
        // case where a projective track from IP reaches the bottom/top faces of the DCH cylinder
        if ((pZ <= m_dchOuterZ / m_dchInnerR * pT) && (pZ >= m_dchOuterZ / m_dchOuterR * pT))
        {
            const float innerExpectedHitRadius(std::max(m_dchInnerR, rInnermostHit));
            const float frac((m_dchOuterZ * pT / pZ - innerExpectedHitRadius) / (m_dchOuterR - innerExpectedHitRadius));
            nExpectedDchHits = frac * m_dchNLayers;
        }

        // calculate number of DCH hits, compare to min number of DCH hits requested
        const int nDchHits(this->GetNDchHits(pTrack));
        const int minDchHits = static_cast<int>(nExpectedDchHits * m_settings.m_minBarrelTrackerHitFractionOfExpected);
        if (nDchHits < minDchHits)
        {
            streamlog_out(WARNING) << " Dropping track : " << momentumAtDca.GetMagnitude() << " Number of DCH hits = " << nDchHits
                                << " < " << minDchHits << std::endl;

            streamlog_out(DEBUG5)  << " track : " << *pTrack
                                << std::endl;
            return false;
        }

        // GM FIXME: make requirement on number of Si wrapper hits configurable
        const int nSiWrapperHits(this->GetNSiWrapperHits(pTrack));
        const int minSiWrapperHits(1);
        if (nSiWrapperHits < minSiWrapperHits) {
            streamlog_out(WARNING) << " Dropping track : " << momentumAtDca.GetMagnitude() << " Number of wrapper hits = " << nDchHits
                                << " < " << minSiWrapperHits << std::endl;
        }
    }
    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

int DDTrackCreatorALLEGRO::GetNVertexHits(const EVENT::Track *const pTrack) const
{
    return pTrack->getSubdetectorHitNumbers()[0] + pTrack->getSubdetectorHitNumbers()[1];  // 0: VXD_barrel; 1: VXD_endcap; 2: DCH; 3: SwiWr_barrel; 4: SiWr_endcap
}

//------------------------------------------------------------------------------------------------------------------------------------------

int DDTrackCreatorALLEGRO::GetNDchHits(const EVENT::Track *const pTrack) const
{
    return pTrack->getSubdetectorHitNumbers()[2];  // 0: VXD_barrel; 1: VXD_endcap; 2: DCH; 3: SwiWr_barrel; 4: SiWr_endcap
}

//------------------------------------------------------------------------------------------------------------------------------------------

int DDTrackCreatorALLEGRO::GetNSiWrapperHits(const EVENT::Track *const pTrack) const
{
    return pTrack->getSubdetectorHitNumbers()[3] + pTrack->getSubdetectorHitNumbers()[4];  // 0: VXD_barrel; 1: VXD_endcap; 2: DCH; 3: SwiWr_barrel; 4: SiWr_endcap
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DDTrackCreatorALLEGRO::DefineTrackPfoUsage(const EVENT::Track *const pTrack, PandoraApi::Track::Parameters &trackParameters) const
{
    bool canFormPfo(false);
    bool canFormClusterlessPfo(false);

    if (trackParameters.m_reachesCalorimeter.Get() && !this->IsParent(pTrack))
    {
        const float d0(std::fabs(pTrack->getD0())), z0(std::fabs(pTrack->getZ0()));
        float rInner(std::numeric_limits<float>::max()), zMin(std::numeric_limits<float>::max());

        // GM: why looping over track hits rather than just taking track state at first hit?
        // this will fail on our tracks currently created from MCParticles without attaching hits to them
        /*
        EVENT::TrackerHitVec trackerHitvec(pTrack->getTrackerHits());
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
        */
        // replace with simpler calculation of rInner and zMin from first hit (though not necessarily the best for e.g. curling track)
        pandora::CartesianVector posAtStart(trackParameters.m_trackStateAtStart.Get().GetPosition());
        rInner = std::sqrt(posAtStart.GetX()*posAtStart.GetX() + posAtStart.GetY()*posAtStart.GetY());
        zMin = std::fabs(posAtStarg.GetZ())

        if (this->PassesQualityCuts(pTrack, trackParameters))
        {
            const pandora::CartesianVector &momentumAtDca(trackParameters.m_momentumAtDca.Get());
            const float pX(momentumAtDca.GetX()), pY(momentumAtDca.GetY()), pZ(momentumAtDca.GetZ());
            const float pT(std::sqrt(pX * pX + pY * pY));

            const float zCutForNonVertexTracks(m_dchInnerR * std::fabs(pZ / pT) + m_settings.m_zCutForNonVertexTracks);
            const bool passRzQualityCuts((zMin < zCutForNonVertexTracks) && (rInner < m_dchInnerR + m_settings.m_maxBarrelTrackerInnerRDistance));

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
                    (rInner < m_dchInnerR + m_settings.m_maxBarrelTrackerInnerRDistance))
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

    // GM debug
    streamlog_out(MESSAGE) << "Track: " << *pTrack << std::endl;
    streamlog_out(MESSAGE) << "Can form PFO: " << canFormPfo << std::endl;
    streamlog_out(MESSAGE) << "Can form clusterless PFO: " << canFormClusterlessPfo << std::endl;
    trackParameters.m_canFormPfo = canFormPfo;
    trackParameters.m_canFormClusterlessPfo = canFormClusterlessPfo;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DDTrackCreatorALLEGRO::TrackReachesECAL(const EVENT::Track *const pTrack, PandoraApi::Track::Parameters &trackParameters) const
{
    // GM debug
    std::cout << "In DDTrackCreatorALLEGRO::TrackReachesECAL" << std::endl;

    // we could probably just return true as done in DDTrackCreatorILD since there are quality checks in DefineTrackPfoUsage()

    // Check that track has hits in wrapper
    if (this->GetNSiWrapperHits(pTrack)<1) {
        std::cout << "Track has no hit in wrapper, thus does not reach ECAL" << std::endl;
        trackParameters.m_reachesCalorimeter = false;
        return;
    }

    // Check that track extrapolates to calo
    pandora::CartesianVector posAtCalo(trackParameters.m_trackStateAtCalorimeter.Get().GetPosition());
    double rAtCalo = std::sqrt(posAtCalo.GetX()*posAtCalo.GetX() + posAtCalo.GetY()*posAtCalo.GetY());
    double zAtCalo = std::fabs(posAtCalo.GetZ());

    // GM debug
    std::cout << "track r, z at ECAL: " << rAtCalo << " , " << zAtCalo << std::endl;

    if (zAtCalo >= m_settings.m_eCalEndCapInnerZ) {
        if (rAtCalo < m_settings.m_eCalEndCapInnerR || rAtCalo > m_settings.m_eCalEndCapOuterR) {
            std::cout << "Track does not reach ECAL endcap" << std::endl;
            trackParameters.m_reachesCalorimeter = false;
            return;
        }
        else {
            std::cout << "Track reaches ECAL endcap" << std::endl;
            trackParameters.m_reachesCalorimeter = true;
            return;
        }
    }
    else {
        if (rAtCalo >= m_settings.m_eCalBarrelInnerR) {
            if (zAtCalo <= m_settings.m_eCalBarrelOuterZ) {
                std::cout << "Track reaches ECAL barrel" << std::endl;
                trackParameters.m_reachesCalorimeter = true;
                return;
            }
        }
        std::cout << "Track does not reach ECAL barrel" << std::endl;
        trackParameters.m_reachesCalorimeter = false;
        return;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DDTrackCreatorALLEGRO::GetTrackStatesAtCalo( EVENT::Track *track,
                                                  lc_content::LCTrackParameters& trackParameters ){
  std::cout << "In DDTrackCreatorALLEGRO::GetTrackStatesAtCalo" << std::endl;
  if( not trackParameters.m_reachesCalorimeter.Get() ) {
      std::cout << "Track does not reach the ECal" << std::endl;
      streamlog_out(DEBUG5) << "Track does not reach the ECal" <<std::endl;
    return;
  }

  const TrackState *trackAtCalo = track->getTrackState(TrackState::AtCalorimeter);
  if( not trackAtCalo ) {
      std::cout << "Track does not have trackState at calorimeter" << std::endl;
      streamlog_out(DEBUG5) << "Track does not have a trackState at calorimeter" <<std::endl;
      streamlog_out(DEBUG3) << toString(track) << std::endl;
      return;
  }

  streamlog_out(DEBUG3) << "Original" << toString(trackAtCalo) << std::endl;

  const auto* tsPosition = trackAtCalo->getReferencePoint();

  if( std::fabs(tsPosition[2]) <  getTrackingRegionExtent()[2] ) {
      std::cout << "Original trackState is at Barrel, should check for extrapolation to Endcap too but won't do" << std::endl;
      streamlog_out(DEBUG5) << "Original trackState is at Barrel" << std::endl;
      pandora::InputTrackState pandoraTrackState;
      this->CopyTrackState( trackAtCalo, pandoraTrackState );
      trackParameters.m_trackStates.push_back( pandoraTrackState );
      // GM FIXME: for the moment, do not extrapolate also to endcap, and return
      return;
  } else { // if track state is in endcap we do not repeat track state calculation, because the barrel cannot be hit
      std::cout << "Original trackState is at Endcap" << std::endl;
      streamlog_out(DEBUG5) << "Original track state is at Endcap" << std::endl;
      pandora::InputTrackState pandoraTrackState;
      this->CopyTrackState( trackAtCalo, pandoraTrackState );
      trackParameters.m_trackStates.push_back( pandoraTrackState );
      return;
  }

  // GM FIXME: for the moment, do not extrapolate also to endcap, and return
  // Wouldn't it be simpler to store second intersection with ECAL (if any, for tracks at ~edge of barrel)
  // to save both track states when doing the original fit (or in TracksFromGenParticlesWithExtrap)?

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
