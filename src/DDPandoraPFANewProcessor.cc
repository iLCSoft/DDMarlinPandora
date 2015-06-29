/**
 *  @file   MarlinPandora/src/DDPandoraPFANewProcessor.cc
 * 
 *  @brief  Implementation of the pandora pfa new processor class.
 * 
 *  $Log: $
 */

#include "marlin/Global.h"
#include "marlin/Exceptions.h"


#include "Api/PandoraApi.h"

#include "LCContent.h"
#if __cplusplus > 199711L
    #include "LCContentFast.h"
#endif

#include "ExternalClusteringAlgorithm.h"
#include "DDPandoraPFANewProcessor.h"

#include <cstdlib>

DDPandoraPFANewProcessor aDDPandoraPFANewProcessor;

DDPandoraPFANewProcessor::PandoraToLCEventMap DDPandoraPFANewProcessor::m_pandoraToLCEventMap;

//------------------------------------------------------------------------------------------------------------------------------------------

DDPandoraPFANewProcessor::DDPandoraPFANewProcessor() :
    Processor("DDPandoraPFANewProcessor"),
    m_pPandora(NULL),
    m_pCaloHitCreator(NULL),
    m_pGeometryCreator(NULL),
    m_pTrackCreator(NULL),
    m_pMCParticleCreator(NULL),
    m_pPfoCreator(NULL)
{
    _description = "Pandora reconstructs clusters and particle flow objects";
    this->ProcessSteeringFile();
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DDPandoraPFANewProcessor::init()
{
    printParameters();
    try
    {
        streamlog_out(MESSAGE) << "DDPandoraPFANewProcessor - Init" << std::endl;
        this->FinaliseSteeringParameters();

        m_pPandora = new pandora::Pandora();
        m_pGeometryCreator = new DDGeometryCreator(m_geometryCreatorSettings, m_pPandora);
        m_pCaloHitCreator = new DDCaloHitCreator(m_caloHitCreatorSettings, m_pPandora);
        m_pTrackCreator = new TrackCreator(m_trackCreatorSettings, m_pPandora);
        m_pMCParticleCreator = new MCParticleCreator(m_mcParticleCreatorSettings, m_pPandora);
        m_pPfoCreator = new PfoCreator(m_pfoCreatorSettings, m_pPandora);

        PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, this->RegisterUserComponents());
        PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, m_pGeometryCreator->CreateGeometry());
        PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::ReadSettings(*m_pPandora, m_settings.m_pandoraSettingsXmlFile));
    }
    catch (pandora::StatusCodeException &statusCodeException)
    {
        streamlog_out(ERROR) << "Failed to initialize marlin pandora: " << statusCodeException.ToString() << std::endl;
        throw statusCodeException;
    }
    catch (std::exception &exception)
    {
        streamlog_out(ERROR) << "Failed to initialize marlin pandora: std exception " << exception.what() << std::endl;
        throw exception;
    }
    catch (...)
    {
        streamlog_out(ERROR) << "Failed to initialize marlin pandora: unrecognized exception" << std::endl;
        throw;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DDPandoraPFANewProcessor::processRunHeader(LCRunHeader *pLCRunHeader)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DDPandoraPFANewProcessor::processEvent(LCEvent *pLCEvent)
{
    try
    {
        streamlog_out(DEBUG) << "DDPandoraPFANewProcessor - Run " << std::endl;
        (void) m_pandoraToLCEventMap.insert(PandoraToLCEventMap::value_type(m_pPandora, pLCEvent));

        PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, m_pMCParticleCreator->CreateMCParticles(pLCEvent));
        PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, m_pTrackCreator->CreateTrackAssociations(pLCEvent));
        PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, m_pTrackCreator->CreateTracks(pLCEvent));
        PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, m_pMCParticleCreator->CreateTrackToMCParticleRelationships(pLCEvent, m_pTrackCreator->GetTrackVector()));
        PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, m_pCaloHitCreator->CreateCaloHits(pLCEvent));
        PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, m_pMCParticleCreator->CreateCaloHitToMCParticleRelationships(pLCEvent, m_pCaloHitCreator->GetCalorimeterHitVector()));

        PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::ProcessEvent(*m_pPandora));
        PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, m_pPfoCreator->CreateParticleFlowObjects(pLCEvent));

        PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::Reset(*m_pPandora));
        this->Reset();
    }
    catch (pandora::StatusCodeException &statusCodeException)
    {
        streamlog_out(ERROR) << "Marlin pandora failed to process event: " << statusCodeException.ToString() << std::endl;
        throw statusCodeException;
    }
    catch (std::exception &exception)
    {
        streamlog_out(ERROR) << "Marlin pandora failed to process event: std exception " << exception.what() << std::endl;
        throw exception;
    }
    catch (EVENT::Exception &exception)
    {
        streamlog_out(ERROR) << "Marlin pandora failed to process event: lcio exception " << exception.what() << std::endl;
        throw exception;
    }
    catch (...)
    {
        streamlog_out(ERROR) << "Marlin pandora failed to process event: unrecognized exception" << std::endl;
        throw;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DDPandoraPFANewProcessor::check(LCEvent */*pLCEvent*/)
{
    streamlog_out(DEBUG) << "DDPandoraPFANewProcessor - Check" << std::endl;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DDPandoraPFANewProcessor::end()
{
    delete m_pPandora;
    delete m_pGeometryCreator;
    delete m_pCaloHitCreator;
    delete m_pTrackCreator;
    delete m_pMCParticleCreator;
    delete m_pPfoCreator;

    streamlog_out(MESSAGE) << "DDPandoraPFANewProcessor - End" << std::endl;
}

//------------------------------------------------------------------------------------------------------------------------------------------

const pandora::Pandora *DDPandoraPFANewProcessor::GetPandora() const
{
    if (NULL == m_pPandora)
        throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_INITIALIZED);

    return m_pPandora;
}

//------------------------------------------------------------------------------------------------------------------------------------------

const EVENT::LCEvent *DDPandoraPFANewProcessor::GetCurrentEvent(const pandora::Pandora *const pPandora)
{
    PandoraToLCEventMap::iterator iter = m_pandoraToLCEventMap.find(pPandora);

    if (m_pandoraToLCEventMap.end() == iter)
        throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_FOUND);

    return iter->second;
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode DDPandoraPFANewProcessor::RegisterUserComponents() const
{
    PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, LCContent::RegisterAlgorithms(*m_pPandora));
#if __cplusplus > 199711L
    PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, LCContentFast::RegisterAlgorithms(*m_pPandora));
#endif

    PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, LCContent::RegisterBasicPlugins(*m_pPandora));

    PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, LCContent::RegisterBFieldPlugin(*m_pPandora,
        m_settings.m_innerBField, m_settings.m_muonBarrelBField, m_settings.m_muonEndCapBField));

    PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, LCContent::RegisterNonLinearityEnergyCorrection(*m_pPandora,
        "NonLinearity", pandora::HADRONIC, m_settings.m_inputEnergyCorrectionPoints, m_settings.m_outputEnergyCorrectionPoints));

    PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::RegisterAlgorithmFactory(*m_pPandora,
        "ExternalClustering", new ExternalClusteringAlgorithm::Factory));

    return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

void DDPandoraPFANewProcessor::ProcessSteeringFile()
{
    registerProcessorParameter("PandoraSettingsXmlFile",
                            "The pandora settings xml file",
                            m_settings.m_pandoraSettingsXmlFile,
                            std::string());

    // Input collections
    registerInputCollections(LCIO::TRACK,
                            "TrackCollections", 
                            "Names of the Track collections used for clustering",
                            m_trackCreatorSettings.m_trackCollections,
                            StringVector());

    registerInputCollections(LCIO::VERTEX,
                            "KinkVertexCollections", 
                            "Name of external kink Vertex collections",
                            m_trackCreatorSettings.m_kinkVertexCollections,
                            StringVector());

    registerInputCollections(LCIO::VERTEX,
                            "ProngVertexCollections", 
                            "Name of external prong Vertex collections",
                            m_trackCreatorSettings.m_prongVertexCollections,
                            StringVector());

    registerInputCollections(LCIO::VERTEX,
                            "SplitVertexCollections", 
                            "Name of external split Vertex collections",
                            m_trackCreatorSettings.m_splitVertexCollections,
                            StringVector());

    registerInputCollections(LCIO::VERTEX,
                            "V0VertexCollections", 
                            "Name of external V0 Vertex collections",
                            m_trackCreatorSettings.m_v0VertexCollections,
                            StringVector());

    registerInputCollections(LCIO::CALORIMETERHIT,
                            "ECalCaloHitCollections", 
                            "Name of the ECAL calo hit collections",
                            m_caloHitCreatorSettings.m_eCalCaloHitCollections,
                            StringVector());

    registerInputCollections(LCIO::CALORIMETERHIT,
                            "HCalCaloHitCollections", 
                            "Name of the HCAL calo hit collections",
                            m_caloHitCreatorSettings.m_hCalCaloHitCollections,
                            StringVector());

    registerInputCollections(LCIO::CALORIMETERHIT,
                            "LCalCaloHitCollections", 
                            "Name of the LCAL calo hit collections",
                            m_caloHitCreatorSettings.m_lCalCaloHitCollections,
                            StringVector());

    registerInputCollections(LCIO::CALORIMETERHIT,
                            "LHCalCaloHitCollections", 
                            "Name of the LHCAL calo hit collections",
                            m_caloHitCreatorSettings.m_lHCalCaloHitCollections,
                            StringVector());

    registerInputCollections(LCIO::CALORIMETERHIT,
                            "MuonCaloHitCollections", 
                            "Name of the muon calo hit collections",
                            m_caloHitCreatorSettings.m_muonCaloHitCollections,
                            StringVector());

    registerInputCollections(LCIO::MCPARTICLE,
                            "MCParticleCollections", 
                            "Name of mc particle collections",
                            m_mcParticleCreatorSettings.m_mcParticleCollections,
                            StringVector());

    registerInputCollections(LCIO::LCRELATION, 
                            "RelCaloHitCollections",
                            "SimCaloHit to CaloHit Relations Collection Name",
                            m_mcParticleCreatorSettings.m_lcCaloHitRelationCollections,
                            StringVector());

    registerInputCollections(LCIO::LCRELATION, 
                            "RelTrackCollections",
                            "Track to MCParticle Relations Collection Name",
                            m_mcParticleCreatorSettings.m_lcTrackRelationCollections,
                            StringVector());

    // Absorber properties
    registerProcessorParameter("AbsorberRadLengthECal",
                            "The absorber radation length in the ECal",
                            m_geometryCreatorSettings.m_absorberRadLengthECal,
                            float(0.2854)); // Default: W, 1 / X0[mm]

    registerProcessorParameter("AbsorberIntLengthECal",
                            "The absorber interaction length in the ECal",
                            m_geometryCreatorSettings.m_absorberIntLengthECal,
                            float(0.0101)); // Default: W, 1 / lambdaI[mm]

    registerProcessorParameter("AbsorberRadLengthHCal",
                            "The absorber radation length in the HCal",
                            m_geometryCreatorSettings.m_absorberRadLengthHCal,
                            float(0.0569)); // Default: Fe, 1 / X0[mm]

    registerProcessorParameter("AbsorberIntLengthHCal",
                            "The absorber interaction length in the HCal",
                            m_geometryCreatorSettings.m_absorberIntLengthHCal,
                            float(0.0060)); // Default: Fe, 1 / lambdaI[mm]

    registerProcessorParameter("AbsorberRadLengthOther",
                            "The absorber radation length in other detector regions",
                            m_geometryCreatorSettings.m_absorberRadLengthOther,
                            float(0.0569)); // Default: Fe, 1 / X0[mm]

    registerProcessorParameter("AbsorberIntLengthOther",
                            "The absorber interaction length in other detector regions",
                            m_geometryCreatorSettings.m_absorberIntLengthOther,
                            float(0.0060)); // Default: Fe, 1 / lambdaI[mm]

    // Name of PFO collection written by MarlinPandora
    registerOutputCollection(LCIO::CLUSTER,
                             "ClusterCollectionName",
                             "Cluster Collection Name",
                             m_pfoCreatorSettings.m_clusterCollectionName,
                             std::string("PandoraPFANewClusters"));

    registerOutputCollection(LCIO::RECONSTRUCTEDPARTICLE,
                             "PFOCollectionName",
                             "PFO Collection Name",
                             m_pfoCreatorSettings.m_pfoCollectionName,
                             std::string("PandoraPFANewPFOs"));

    // Calibration constants
    registerProcessorParameter("ECalToMipCalibration",
                            "The calibration from deposited ECal energy to mip",
                            m_caloHitCreatorSettings.m_eCalToMip,
                            float(1.));

    registerProcessorParameter("HCalToMipCalibration",
                            "The calibration from deposited HCal energy to mip",
                            m_caloHitCreatorSettings.m_hCalToMip,
                            float(1.));

    registerProcessorParameter("ECalMipThreshold",
                            "Threshold for creating calo hits in the ECal, units mip",
                            m_caloHitCreatorSettings.m_eCalMipThreshold,
                            float(0.));

    registerProcessorParameter("MuonToMipCalibration",
                            "The calibration from deposited Muon energy to mip",
                            m_caloHitCreatorSettings.m_muonToMip,
                            float(1.));

    registerProcessorParameter("HCalMipThreshold",
                            "Threshold for creating calo hits in the HCal, units mip",
                            m_caloHitCreatorSettings.m_hCalMipThreshold,
                            float(0.));

    registerProcessorParameter("ECalToEMGeVCalibration",
                            "The calibration from deposited ECal energy to EM energy",
                            m_caloHitCreatorSettings.m_eCalToEMGeV,
                            float(1.));

    registerProcessorParameter("HCalToEMGeVCalibration",
                            "The calibration from deposited HCal energy to EM energy",
                            m_caloHitCreatorSettings.m_hCalToEMGeV,
                            float(1.));

    registerProcessorParameter("ECalToHadGeVCalibrationEndCap",
                            "The calibration from deposited ECal energy to hadronic energy",
                            m_caloHitCreatorSettings.m_eCalToHadGeVEndCap,
                            float(1.));

    registerProcessorParameter("ECalToHadGeVCalibrationBarrel",
                            "The calibration from deposited ECal energy to hadronic energy",
                            m_caloHitCreatorSettings.m_eCalToHadGeVBarrel,
                            float(1.));

    registerProcessorParameter("HCalToHadGeVCalibration",
                            "The calibration from deposited HCal energy to hadronic energy",
                            m_caloHitCreatorSettings.m_hCalToHadGeV,
                            float(1.));

    registerProcessorParameter("DigitalMuonHits",
                            "Treat muon hits as digital",
                            m_caloHitCreatorSettings.m_muonDigitalHits,
                            int(1));

    registerProcessorParameter("MuonHitEnergy",
                            "The energy for a digital muon calorimeter hit, units GeV",
                            m_caloHitCreatorSettings.m_muonHitEnergy,
                            float(0.5));

    registerProcessorParameter("MaxHCalHitHadronicEnergy",
                            "The maximum hadronic energy allowed for a single hcal hit",
                            m_caloHitCreatorSettings.m_maxHCalHitHadronicEnergy,
                            float(10000.));

    registerProcessorParameter("NOuterSamplingLayers",
                            "Number of layers from edge for hit to be flagged as an outer layer hit",
                            m_caloHitCreatorSettings.m_nOuterSamplingLayers,
                            int(3));

    registerProcessorParameter("LayersFromEdgeMaxRearDistance",
                            "Maximum number of layers from candidate outer layer hit to rear of detector",
                            m_caloHitCreatorSettings.m_layersFromEdgeMaxRearDistance,
                            float(250.f));

    // B-field parameters
    registerProcessorParameter("MuonBarrelBField",
                            "The bfield in the muon barrel, units Tesla",
                            m_settings.m_muonBarrelBField,
                            float(-1.5f));

    registerProcessorParameter("MuonEndCapBField",
                            "The bfield in the muon endcap, units Tesla",
                            m_settings.m_muonEndCapBField,
                            float(0.01f));

    // Track relationship parameters
    registerProcessorParameter("ShouldFormTrackRelationships",
                            "Whether to form pandora track relationships using v0 and kink info",
                            m_trackCreatorSettings.m_shouldFormTrackRelationships,
                            int(1));

    // Initial track hit specifications
   registerProcessorParameter("MinTrackHits",
                            "Track quality cut: the minimum number of track hits",
                            m_trackCreatorSettings.m_minTrackHits,
                            int(5));

   registerProcessorParameter("MinFtdTrackHits",
                            "Track quality cut: the minimum number of ftd track hits for ftd only tracks",
                            m_trackCreatorSettings.m_minFtdTrackHits,
                            int(0));

   registerProcessorParameter("MaxTrackHits",
                            "Track quality cut: the maximum number of track hits",
                            m_trackCreatorSettings.m_maxTrackHits,
                            int(5000));

   registerProcessorParameter("UseOldTrackStateCalculation",
                            "Whether to calculate track states manually, rather than copy stored fitter values",
                            m_trackCreatorSettings.m_useOldTrackStateCalculation,
                            int(0));

    // Track PFO usage parameters
    registerProcessorParameter("D0TrackCut",
                            "Track d0 cut used to determine whether track can be used to form pfo",
                            m_trackCreatorSettings.m_d0TrackCut,
                            float(50.));

    registerProcessorParameter("Z0TrackCut",
                            "Track z0 cut used to determine whether track can be used to form pfo",
                            m_trackCreatorSettings.m_z0TrackCut,
                            float(50.));

    registerProcessorParameter("UseNonVertexTracks",
                            "Whether can form pfos from tracks that don't start at vertex",
                            m_trackCreatorSettings.m_usingNonVertexTracks,
                            int(1));

    registerProcessorParameter("UseUnmatchedNonVertexTracks",
                            "Whether can form pfos from unmatched tracks that don't start at vertex",
                            m_trackCreatorSettings.m_usingUnmatchedNonVertexTracks,
                            int(0));

    registerProcessorParameter("UseUnmatchedVertexTracks",
                            "Whether can form pfos from unmatched tracks that start at vertex",
                            m_trackCreatorSettings.m_usingUnmatchedVertexTracks,
                            int(1));

    registerProcessorParameter("UnmatchedVertexTrackMaxEnergy",
                            "Maximum energy for unmatched vertex track",
                            m_trackCreatorSettings.m_unmatchedVertexTrackMaxEnergy,
                            float(5.));

    registerProcessorParameter("D0UnmatchedVertexTrackCut",
                            "d0 cut used to determine whether unmatched vertex track can form pfo",
                            m_trackCreatorSettings.m_d0UnmatchedVertexTrackCut,
                            float(5.));

    registerProcessorParameter("Z0UnmatchedVertexTrackCut",
                            "z0 cut used to determine whether unmatched vertex track can form pfo",
                            m_trackCreatorSettings.m_z0UnmatchedVertexTrackCut,
                            float(5.));

    registerProcessorParameter("ZCutForNonVertexTracks",
                            "Non vtx track z cut to determine whether track can be used to form pfo",
                            m_trackCreatorSettings.m_zCutForNonVertexTracks,
                            float(250.));

    // Track "reaches ecal" parameters
    registerProcessorParameter("ReachesECalNTpcHits",
                            "Minimum number of tpc hits to consider track as reaching ecal",
                            m_trackCreatorSettings.m_reachesECalNTpcHits,
                            int(11));

    registerProcessorParameter("ReachesECalNFtdHits",
                            "Minimum number of ftd hits to consider track as reaching ecal",
                            m_trackCreatorSettings.m_reachesECalNFtdHits,
                            int(4));

    registerProcessorParameter("ReachesECalTpcOuterDistance",
                            "Max distance from track to tpc r max to id whether track reaches ecal",
                            m_trackCreatorSettings.m_reachesECalTpcOuterDistance,
                            float(-100.));

    registerProcessorParameter("ReachesECalMinFtdLayer",
                            "Min FTD layer for track to be considered to have reached ecal",
                            m_trackCreatorSettings.m_reachesECalMinFtdLayer,
                            int(9));

    registerProcessorParameter("ReachesECalTpcZMaxDistance",
                            "Max distance from track to tpc z max to id whether track reaches ecal",
                            m_trackCreatorSettings.m_reachesECalTpcZMaxDistance,
                            float(-50.));

    registerProcessorParameter("ReachesECalFtdZMaxDistance",
                            "Max distance from track hit to ftd z position to identify ftd hits",
                            m_trackCreatorSettings.m_reachesECalFtdZMaxDistance,
                            float(1.));

    registerProcessorParameter("CurvatureToMomentumFactor",
                            "Constant relating track curvature in b field to momentum",
                            m_trackCreatorSettings.m_curvatureToMomentumFactor,
                            float(0.3 / 2000.));

    registerProcessorParameter("MinTrackECalDistanceFromIp",
                            "Sanity check on separation between ip and track projected ecal position",
                            m_trackCreatorSettings.m_minTrackECalDistanceFromIp,
                            float(100.));

    // Final track quality parameters
    registerProcessorParameter("MaxTrackSigmaPOverP",
                            "Cut on fractional track momentum error",
                            m_trackCreatorSettings.m_maxTrackSigmaPOverP,
                            float(0.15));

    registerProcessorParameter("MinMomentumForTrackHitChecks",
                            "Min track momentum required to perform final quality checks on number of hits",
                            m_trackCreatorSettings.m_minMomentumForTrackHitChecks,
                            float(1.));

    registerProcessorParameter("TpcMembraneMaxZ",
                            "Tpc membrane max z coordinate",
                            m_trackCreatorSettings.m_tpcMembraneMaxZ,
                            float(10.));

    registerProcessorParameter("MinTpcHitFractionOfExpected",
                            "Cut on fractional of expected number of TPC hits",
                            m_trackCreatorSettings.m_minTpcHitFractionOfExpected,
                            float(0.20));

    registerProcessorParameter("MinFtdHitsForTpcHitFraction",
                            "Cut on minimum number of FTD hits of TPC hit fraction to be applied",
                            m_trackCreatorSettings.m_minFtdHitsForTpcHitFraction,
                            int(2));

    registerProcessorParameter("MaxTpcInnerRDistance",
                            "Track cut on distance from tpc inner r to id whether track can form pfo",
                            m_trackCreatorSettings.m_maxTpcInnerRDistance,
                            float(50.));

    // Additional geometry parameters
    registerProcessorParameter("ECalEndCapInnerSymmetryOrder",
                            "ECal end cap inner symmetry order (missing from ILD gear files)",
                            m_geometryCreatorSettings.m_eCalEndCapInnerSymmetryOrder,
                            int(4));

    registerProcessorParameter("ECalEndCapInnerPhiCoordinate",
                            "ECal end cap inner phi coordinate (missing from ILD gear files)",
                            m_geometryCreatorSettings.m_eCalEndCapInnerPhiCoordinate,
                            float(0.));

    registerProcessorParameter("ECalEndCapOuterSymmetryOrder",
                            "ECal end cap outer symmetry order (missing from ILD gear files)",
                            m_geometryCreatorSettings.m_eCalEndCapOuterSymmetryOrder,
                            int(8));

    registerProcessorParameter("ECalEndCapOuterPhiCoordinate",
                            "ECal end cap outer phi coordinate (missing from ILD gear files)",
                            m_geometryCreatorSettings.m_eCalEndCapOuterPhiCoordinate,
                            float(0.));

    registerProcessorParameter("HCalEndCapInnerSymmetryOrder",
                            "HCal end cap inner symmetry order (missing from ILD gear files)",
                            m_geometryCreatorSettings.m_hCalEndCapInnerSymmetryOrder,
                            int(4));

    registerProcessorParameter("HCalEndCapInnerPhiCoordinate",
                            "HCal end cap inner phi coordinate (missing from ILD gear files)",
                            m_geometryCreatorSettings.m_hCalEndCapInnerPhiCoordinate,
                            float(0.));

    registerProcessorParameter("HCalEndCapOuterSymmetryOrder",
                            "HCal end cap outer symmetry order (missing from ILD gear files)",
                            m_geometryCreatorSettings.m_hCalEndCapOuterSymmetryOrder,
                            int(16));

    registerProcessorParameter("HCalEndCapOuterPhiCoordinate",
                            "HCal end cap outer phi coordinate (missing from ILD gear files)",
                            m_geometryCreatorSettings.m_hCalEndCapOuterPhiCoordinate,
                            float(0.));

    registerProcessorParameter("HCalRingInnerSymmetryOrder",
                            "HCal ring inner symmetry order (missing from ILD gear files)",
                            m_geometryCreatorSettings.m_hCalRingInnerSymmetryOrder,
                            int(8));

    registerProcessorParameter("HCalRingInnerPhiCoordinate",
                            "HCal ring inner phi coordinate (missing from ILD gear files)",
                            m_geometryCreatorSettings.m_hCalRingInnerPhiCoordinate,
                            float(0.));

    registerProcessorParameter("HCalRingOuterSymmetryOrder",
                            "HCal ring outer symmetry order (missing from ILD gear files)",
                            m_geometryCreatorSettings.m_hCalRingOuterSymmetryOrder,
                            int(16));

    registerProcessorParameter("HCalRingOuterPhiCoordinate",
                            "HCal ring outer phi coordinate (missing from ILD gear files)",
                            m_geometryCreatorSettings.m_hCalRingOuterPhiCoordinate,
                            float(0.));

    // For Strip Splitting method and also for hybrid ECAL
    registerProcessorParameter("StripSplittingOn",
                            "To use strip splitting algorithm, this should be true",
                            m_caloHitCreatorSettings.m_stripSplittingOn,
                            int(0));

    // For Strip Splitting method and also for hybrid ECAL
    registerProcessorParameter("UseEcalScLayers",
                            "To use scintillator layers ~ hybrid ECAL, this should be true",
                            m_caloHitCreatorSettings.m_useEcalScLayers,
                            int(0));

    // Parameters for hybrid ECAL
    // Energy to MIP for Si-layers and Sc-layers, respectively.
    //Si
    registerProcessorParameter("ECalSiToMipCalibration",
                            "The calibration from deposited Si-layer energy to mip",
                            m_caloHitCreatorSettings.m_eCalSiToMip,
                            float(1.));

    //Sc
    registerProcessorParameter("ECalScToMipCalibration",
                            "The calibration from deposited Sc-layer energy to mip",
                            m_caloHitCreatorSettings.m_eCalScToMip,
                            float(1.));

    // MipThreshold for Si-layers and Sc-layers, respectively.
    // Si
    registerProcessorParameter("ECalSiMipThreshold",
                            "Threshold for creating calo hits in the Si-layers of ECAL, units mip",
                            m_caloHitCreatorSettings.m_eCalSiMipThreshold,
                            float(0.));

    //Sc
    registerProcessorParameter("ECalScMipThreshold",
                            "Threshold for creating calo hits in the Sc-layers of ECAL, units mip",
                            m_caloHitCreatorSettings.m_eCalScMipThreshold,
                            float(0.));

    // EcalToEM for Si-layers and Sc-layers, respectively.
    //Si
    registerProcessorParameter("ECalSiToEMGeVCalibration",
                            "The calibration from deposited Si-layer energy to EM energy",
                            m_caloHitCreatorSettings.m_eCalSiToEMGeV,
                            float(1.));

    //Sc
    registerProcessorParameter("ECalScToEMGeVCalibration",
                            "The calibration from deposited Sc-layer energy to EM energy",
                            m_caloHitCreatorSettings.m_eCalScToEMGeV,
                            float(1.));

    // EcalToHad for Si-layers and Sc-layers of the endcaps, respectively.
    //Si
    registerProcessorParameter("ECalSiToHadGeVCalibrationEndCap",
                            "The calibration from deposited Si-layer energy on the enecaps to hadronic energy",
                            m_caloHitCreatorSettings.m_eCalSiToHadGeVEndCap,
                            float(1.));

    //Sc
    registerProcessorParameter("ECalScToHadGeVCalibrationEndCap",
                            "The calibration from deposited Sc-layer energy on the endcaps to hadronic energy",
                            m_caloHitCreatorSettings.m_eCalScToHadGeVEndCap,
                            float(1.));

    // EcalToHad for Si-layers and Sc-layers of the barrel, respectively.
    //Si
    registerProcessorParameter("ECalSiToHadGeVCalibrationBarrel",
                            "The calibration from deposited Si-layer energy on the barrel to hadronic energy",
                            m_caloHitCreatorSettings.m_eCalSiToHadGeVBarrel,
                            float(1.));

    //Sc
    registerProcessorParameter("ECalScToHadGeVCalibrationBarrel",
                            "The calibration from deposited Sc-layer energy to the barrel hadronic energy",
                            m_caloHitCreatorSettings.m_eCalScToHadGeVBarrel,
                            float(1.));

    // Hadronic energy non-linearity correction
    registerProcessorParameter("InputEnergyCorrectionPoints",
                            "The input energy points for hadronic energy correction",
                            m_settings.m_inputEnergyCorrectionPoints,
                            FloatVector());

    registerProcessorParameter("OutputEnergyCorrectionPoints",
                            "The output energy points for hadronic energy correction",
                            m_settings.m_outputEnergyCorrectionPoints,
                            FloatVector());
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DDPandoraPFANewProcessor::FinaliseSteeringParameters()
{
    // ATTN: This function seems to be necessary for operations that cannot easily be performed at construction of the processor,
    // when the steering file is parsed e.g. the call to GEAR to get the inner bfield
    m_caloHitCreatorSettings.m_absorberRadLengthECal = m_geometryCreatorSettings.m_absorberRadLengthECal;
    m_caloHitCreatorSettings.m_absorberIntLengthECal = m_geometryCreatorSettings.m_absorberIntLengthECal;
    m_caloHitCreatorSettings.m_absorberRadLengthHCal = m_geometryCreatorSettings.m_absorberRadLengthHCal;
    m_caloHitCreatorSettings.m_absorberIntLengthHCal = m_geometryCreatorSettings.m_absorberIntLengthHCal;
    m_caloHitCreatorSettings.m_absorberRadLengthOther = m_geometryCreatorSettings.m_absorberRadLengthOther;
    m_caloHitCreatorSettings.m_absorberIntLengthOther = m_geometryCreatorSettings.m_absorberIntLengthOther;

    m_caloHitCreatorSettings.m_hCalEndCapInnerSymmetryOrder = m_geometryCreatorSettings.m_hCalEndCapInnerSymmetryOrder;
    m_caloHitCreatorSettings.m_hCalEndCapInnerPhiCoordinate = m_geometryCreatorSettings.m_hCalEndCapInnerPhiCoordinate;

    m_trackCreatorSettings.m_prongSplitVertexCollections = m_trackCreatorSettings.m_prongVertexCollections;
    m_trackCreatorSettings.m_prongSplitVertexCollections.insert(m_trackCreatorSettings.m_prongSplitVertexCollections.end(),
        m_trackCreatorSettings.m_splitVertexCollections.begin(), m_trackCreatorSettings.m_splitVertexCollections.end());

    ///FIXME! CHANGE TO LCDD 

    std::cout<<"BIG WARNING!: HARD CODED BFIELD VALUE TO BE ACCESSED FROM LCDD!"<<std::endl;
    m_settings.m_innerBField = 4.0;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DDPandoraPFANewProcessor::Reset()
{
    m_pCaloHitCreator->Reset();
    m_pTrackCreator->Reset();

    PandoraToLCEventMap::iterator iter = m_pandoraToLCEventMap.find(m_pPandora);

    if (m_pandoraToLCEventMap.end() == iter)
        throw pandora::StatusCodeException(pandora::STATUS_CODE_FAILURE);

    m_pandoraToLCEventMap.erase(iter);
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

DDPandoraPFANewProcessor::Settings::Settings() :
    m_innerBField(3.5f),
    m_muonBarrelBField(-1.5f),
    m_muonEndCapBField(0.01f)
{
}
