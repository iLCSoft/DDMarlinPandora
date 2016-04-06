/**
 *  @file   DDMarlinPandora/src/DDPandoraPFANewProcessor.cc
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

#include "DDExternalClusteringAlgorithm.h"
#include "DDPandoraPFANewProcessor.h"

#include "DD4hep/LCDD.h"
#include "DD4hep/DD4hepUnits.h"
#include "DD4hep/DetType.h"
#include "DDRec/DetectorData.h"
#include "DD4hep/DetectorSelector.h"

#include "DDTrackCreatorILD.h"
#include "DDTrackCreatorCLIC.h"


#include <cstdlib>

DDPandoraPFANewProcessor aDDPandoraPFANewProcessor;

double getFieldFromLCDD(){
  
  DD4hep::Geometry::LCDD& lcdd = DD4hep::Geometry::LCDD::getInstance();
  const double position[3]={0,0,0}; // position to calculate magnetic field at (the origin in this case)
  double magneticFieldVector[3]={0,0,0}; // initialise object to hold magnetic field
  lcdd.field().magneticField(position,magneticFieldVector); // get the magnetic field vector from DD4hep
  
  return magneticFieldVector[2]/dd4hep::tesla; // z component at (0,0,0)
  
}

//Not needed anymore; to be removed
// double getCoilOuterR(){
//     
//   try{
//     DD4hep::Geometry::LCDD & lcdd = DD4hep::Geometry::LCDD::getInstance();
//     const std::vector< DD4hep::Geometry::DetElement>& theDetectors = DD4hep::Geometry::DetectorSelector(lcdd).detectors(  DD4hep::DetType::COIL ) ;
//     //access the detelement and create a shape from the envelope since only minimal info needed
//     DD4hep::Geometry::Tube coilTube = DD4hep::Geometry::Tube( theDetectors.at(0).volume().solid() )  ;
//     return coilTube->GetRmax()/ dd4hep::mm;
//   } catch ( std::exception & e ) {
//       
//           streamlog_out(ERROR)<< "BIG WARNING! CANNOT GET EXTENSION FOR COIL: "<<e.what()<<" MAKE SURE YOU CHANGE THIS!"<< std::endl;
// 
//   }
//   
//   return 0;
// }


///Not needed anymore. To be removed
// DD4hep::DDRec::LayeredCalorimeterData * getExtension(std::string detectorName){
//   
//   
//   DD4hep::DDRec::LayeredCalorimeterData * theExtension = 0;
//   
//   try {
//     DD4hep::Geometry::LCDD & lcdd = DD4hep::Geometry::LCDD::getInstance();
//     const DD4hep::Geometry::DetElement & theDetector = lcdd.detector(detectorName);
//     theExtension = theDetector.extension<DD4hep::DDRec::LayeredCalorimeterData>();
//     //     std::cout<< "DEBUG: in getExtension(\""<<detectorName<<"\"): size of layers: "<<theExtension->layers.size()<<" positions not shown. "<<std::endl;
//     
//     //     for(int i=0; i< theExtension->layers.size(); i++){
//     //       std::cout<<theExtension->layers[i].distance/dd4hep::mm<<" ";
//     //     }
//     //     std::cout<<std::endl;
//   } catch ( ... ){
//     
//     streamlog_out(ERROR) << "BIG WARNING! EXTENSION DOES NOT EXIST FOR " << detectorName<<". MAKE SURE YOU CHANGE THIS!"<< std::endl;
//   }
//   
//   return theExtension;
// }


DD4hep::DDRec::LayeredCalorimeterData * getExtension(unsigned int includeFlag, unsigned int excludeFlag=0) {
  
  
  DD4hep::DDRec::LayeredCalorimeterData * theExtension = 0;
  
  DD4hep::Geometry::LCDD & lcdd = DD4hep::Geometry::LCDD::getInstance();
  const std::vector< DD4hep::Geometry::DetElement>& theDetectors = DD4hep::Geometry::DetectorSelector(lcdd).detectors(  includeFlag, excludeFlag ) ;
  
  
  streamlog_out( DEBUG2 ) << " getExtension :  includeFlag: " << DD4hep::DetType( includeFlag ) << " excludeFlag: " << DD4hep::DetType( excludeFlag ) 
			  << "  found : " << theDetectors.size() << "  - first det: " << theDetectors.at(0).name() << std::endl ;
  
  if( theDetectors.size()  != 1 ){
    
    std::stringstream es ;
    es << " getExtension: selection is not unique (or empty)  includeFlag: " << DD4hep::DetType( includeFlag ) << " excludeFlag: " << DD4hep::DetType( excludeFlag ) 
       << " --- found detectors : " ;
    for( unsigned i=0, N= theDetectors.size(); i<N ; ++i ){
      es << theDetectors.at(i).name() << ", " ; 
    }
    throw std::runtime_error( es.str() ) ;
  }
  
  theExtension = theDetectors.at(0).extension<DD4hep::DDRec::LayeredCalorimeterData>();
  
  return theExtension;
}

std::vector<double> getTrackingRegionExtent(){
  
  ///Rmin, Rmax, Zmax
  std::vector<double> extent;
  
  extent.reserve(3);
  
  DD4hep::Geometry::LCDD & lcdd = DD4hep::Geometry::LCDD::getInstance();
  
  
  
  extent[0]=0.1; ///FIXME! CLIC-specific: Inner radius was set to 0 for SiD-type detectors
  extent[1]=lcdd.constantAsDouble("tracker_region_rmax")/dd4hep::mm;
  extent[2]=lcdd.constantAsDouble("tracker_region_zmax")/dd4hep::mm;

  return extent;
  
  
}

DDPandoraPFANewProcessor::PandoraToLCEventMap DDPandoraPFANewProcessor::m_pandoraToLCEventMap;

//------------------------------------------------------------------------------------------------------------------------------------------

DDPandoraPFANewProcessor::DDPandoraPFANewProcessor() :
    Processor("DDPandoraPFANewProcessor"),
    m_pPandora(NULL),
    m_pCaloHitCreator(NULL),
    m_pGeometryCreator(NULL),
    m_pTrackCreator(NULL),
    m_pDDMCParticleCreator(NULL),
    m_pDDPfoCreator(NULL)
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
        
        
        ///FIXME: IMPLEMENT FACTORY
        if (m_settings.m_trackCreatorName == "DDTrackCreatorCLIC")
            m_pTrackCreator = new DDTrackCreatorCLIC(m_trackCreatorSettings, m_pPandora);
        else if (m_settings.m_trackCreatorName == "DDTrackCreatorILD")
            m_pTrackCreator = new DDTrackCreatorILD(m_trackCreatorSettings, m_pPandora);
        else
            streamlog_out(ERROR) << "Unknown DDTrackCreator: "<<m_settings.m_trackCreatorName << std::endl;


        m_pDDMCParticleCreator = new DDMCParticleCreator(m_mcParticleCreatorSettings, m_pPandora);
        m_pDDPfoCreator = new DDPfoCreator(m_pfoCreatorSettings, m_pPandora);

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

        PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, m_pDDMCParticleCreator->CreateMCParticles(pLCEvent));
        PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, m_pTrackCreator->CreateTrackAssociations(pLCEvent));
        PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, m_pTrackCreator->CreateTracks(pLCEvent));
        PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, m_pDDMCParticleCreator->CreateTrackToMCParticleRelationships(pLCEvent, m_pTrackCreator->GetTrackVector()));
        PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, m_pCaloHitCreator->CreateCaloHits(pLCEvent));
        PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, m_pDDMCParticleCreator->CreateCaloHitToMCParticleRelationships(pLCEvent, m_pCaloHitCreator->GetCalorimeterHitVector()));

        PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::ProcessEvent(*m_pPandora));
        PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, m_pDDPfoCreator->CreateParticleFlowObjects(pLCEvent));

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
    delete m_pDDMCParticleCreator;
    delete m_pDDPfoCreator;

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
        "ExternalClustering", new DDExternalClusteringAlgorithm::Factory));

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

       registerProcessorParameter("CreateGaps",
                            "Decides whether to create gaps in the geometry (ILD-specific)",
                            m_geometryCreatorSettings.m_createGaps,
                            bool(true)); 
    
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
    registerOutputCollection(LCIO::VERTEX,
                             "StartVertexCollectionName",
                             "Start Vertex Collection Name",
                             m_pfoCreatorSettings.m_startVertexCollectionName,
                             std::string("PandoraPFANewStartVertices"));

    registerProcessorParameter("StartVertexAlgorithmName",
                            "The algorithm name for filling start vertex",
                            m_pfoCreatorSettings.m_startVertexAlgName,
                            std::string("PandoraPFANew"));

    // Energy resolution parameters
    registerProcessorParameter("EMStochasticTerm",
                            "The stochastic term for EM shower",
                            m_pfoCreatorSettings.m_emStochasticTerm,
                            float(0.17));

    registerProcessorParameter("HadStochasticTerm",
                            "The stochastic term for Hadronic shower",
                            m_pfoCreatorSettings.m_hadStochasticTerm,
                            float(0.6));

    registerProcessorParameter("EMConstantTerm",
                            "The constant term for EM shower",
                            m_pfoCreatorSettings.m_emConstantTerm,
                            float(0.01));

    registerProcessorParameter("HadConstantTerm",
                            "The constant term for Hadronic shower",
                            m_pfoCreatorSettings.m_hadConstantTerm,
                            float(0.03));

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
    registerProcessorParameter("ReachesECalNBarrelTrackerHits",
                            "Minimum number of BarrelTracker hits to consider track as reaching ecal",
                            m_trackCreatorSettings.m_reachesECalNBarrelTrackerHits,
                            int(11));

    registerProcessorParameter("ReachesECalNFtdHits",
                            "Minimum number of ftd hits to consider track as reaching ecal",
                            m_trackCreatorSettings.m_reachesECalNFtdHits,
                            int(4));

    registerProcessorParameter("ReachesECalBarrelTrackerOuterDistance",
                            "Max distance from track to BarrelTracker r max to id whether track reaches ecal",
                            m_trackCreatorSettings.m_reachesECalBarrelTrackerOuterDistance,
                            float(-100.));

    registerProcessorParameter("ReachesECalMinFtdLayer",
                            "Min FTD layer for track to be considered to have reached ecal",
                            m_trackCreatorSettings.m_reachesECalMinFtdLayer,
                            int(9));

    registerProcessorParameter("ReachesECalBarrelTrackerZMaxDistance",
                            "Max distance from track to BarrelTracker z max to id whether track reaches ecal",
                            m_trackCreatorSettings.m_reachesECalBarrelTrackerZMaxDistance,
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

    registerProcessorParameter("MinBarrelTrackerHitFractionOfExpected",
                            "Cut on fractional of expected number of BarrelTracker hits",
                            m_trackCreatorSettings.m_minBarrelTrackerHitFractionOfExpected,
                            float(0.20));

    registerProcessorParameter("MinFtdHitsForBarrelTrackerHitFraction",
                            "Cut on minimum number of FTD hits of BarrelTracker hit fraction to be applied",
                            m_trackCreatorSettings.m_minFtdHitsForBarrelTrackerHitFraction,
                            int(2));

    registerProcessorParameter("MaxBarrelTrackerInnerRDistance",
                            "Track cut on distance from BarrelTracker inner r to id whether track can form pfo",
                            m_trackCreatorSettings.m_maxBarrelTrackerInnerRDistance,
                            float(50.));

   
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
    
    
    ///EXTRA PARAMETERS FROM NIKIFOROS
    registerProcessorParameter("TrackCreatorName",
                               "The name of the DDTrackCreator implementation",
                               m_settings.m_trackCreatorName,
                               std::string("DDTrackCreatorCLIC")); 
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DDPandoraPFANewProcessor::FinaliseSteeringParameters()
{
    // ATTN: This function seems to be necessary for operations that cannot easily be performed at construction of the processor,
    // when the steering file is parsed e.g. the call to GEAR to get the inner bfield
    m_trackCreatorSettings.m_prongSplitVertexCollections = m_trackCreatorSettings.m_prongVertexCollections;
    m_trackCreatorSettings.m_prongSplitVertexCollections.insert(m_trackCreatorSettings.m_prongSplitVertexCollections.end(),m_trackCreatorSettings.m_splitVertexCollections.begin(),m_trackCreatorSettings.m_splitVertexCollections.end());
    
    m_trackCreatorSettings.m_bField=getFieldFromLCDD();
    
    //Get ECal Barrel extension by type, ignore plugs and rings 
    const DD4hep::DDRec::LayeredCalorimeterData * eCalBarrelExtension= getExtension( ( DD4hep::DetType::CALORIMETER | DD4hep::DetType::ELECTROMAGNETIC | DD4hep::DetType::BARREL), 
										     ( DD4hep::DetType::AUXILIARY  |  DD4hep::DetType::FORWARD ) );
    //Get ECal Endcap extension by type, ignore plugs and rings 
    const DD4hep::DDRec::LayeredCalorimeterData * eCalEndcapExtension= getExtension( ( DD4hep::DetType::CALORIMETER | DD4hep::DetType::ELECTROMAGNETIC | DD4hep::DetType::ENDCAP), 
										     ( DD4hep::DetType::AUXILIARY |  DD4hep::DetType::FORWARD  ) );
    //Get HCal Barrel extension by type, ignore plugs and rings 
    const DD4hep::DDRec::LayeredCalorimeterData * hCalBarrelExtension= getExtension( ( DD4hep::DetType::CALORIMETER | DD4hep::DetType::HADRONIC | DD4hep::DetType::BARREL), 
										     ( DD4hep::DetType::AUXILIARY |  DD4hep::DetType::FORWARD ) );
      //Get HCal Endcap extension by type, ignore plugs and rings 
    const DD4hep::DDRec::LayeredCalorimeterData * hCalEndcapExtension= getExtension( ( DD4hep::DetType::CALORIMETER | DD4hep::DetType::HADRONIC | DD4hep::DetType::ENDCAP), 
										     ( DD4hep::DetType::AUXILIARY |  DD4hep::DetType::FORWARD ) );
    //Get Muon Barrel extension by type, ignore plugs and rings 
    const DD4hep::DDRec::LayeredCalorimeterData * muonBarrelExtension= getExtension( ( DD4hep::DetType::CALORIMETER | DD4hep::DetType::MUON | DD4hep::DetType::BARREL), 
										     ( DD4hep::DetType::AUXILIARY |  DD4hep::DetType::FORWARD ) );
    //fg: muon endcap is not used :
    // //Get Muon Endcap extension by type, ignore plugs and rings 
    // const DD4hep::DDRec::LayeredCalorimeterData * muonEndcapExtension= getExtension( ( DD4hep::DetType::CALORIMETER | DD4hep::DetType::MUON | DD4hep::DetType::ENDCAP), ( DD4hep::DetType::AUXILIARY ) );   
    
    //Get COIL extension
    const DD4hep::DDRec::LayeredCalorimeterData * coilExtension= getExtension( ( DD4hep::DetType::COIL ) );   
  
    
    m_trackCreatorSettings.m_eCalBarrelInnerSymmetry        =   eCalBarrelExtension->inner_symmetry;
    m_trackCreatorSettings.m_eCalBarrelInnerPhi0            =   eCalBarrelExtension->inner_phi0/dd4hep::rad;
    m_trackCreatorSettings.m_eCalBarrelInnerR               =   eCalBarrelExtension->extent[0]/dd4hep::mm;
    m_trackCreatorSettings.m_eCalEndCapInnerZ               =   eCalEndcapExtension->extent[2]/dd4hep::mm;
                                                            
    m_caloHitCreatorSettings.m_eCalBarrelOuterZ             =   eCalBarrelExtension->extent[3]/dd4hep::mm;
    m_caloHitCreatorSettings.m_hCalBarrelOuterZ             =   hCalBarrelExtension->extent[3]/dd4hep::mm;
    m_caloHitCreatorSettings.m_muonBarrelOuterZ             =   muonBarrelExtension->extent[3]/dd4hep::mm;
    m_caloHitCreatorSettings.m_coilOuterR                   =   coilExtension->extent[1]/dd4hep::mm;
    m_caloHitCreatorSettings.m_eCalBarrelInnerPhi0          =   eCalBarrelExtension->inner_phi0/dd4hep::rad;
    m_caloHitCreatorSettings.m_eCalBarrelInnerSymmetry      =   eCalBarrelExtension->inner_symmetry;
    m_caloHitCreatorSettings.m_hCalBarrelInnerPhi0          =   hCalBarrelExtension->inner_phi0/dd4hep::rad;
    m_caloHitCreatorSettings.m_hCalBarrelInnerSymmetry      =   hCalBarrelExtension->inner_symmetry;
    m_caloHitCreatorSettings.m_muonBarrelInnerPhi0          =   muonBarrelExtension->inner_phi0/dd4hep::rad;
    m_caloHitCreatorSettings.m_muonBarrelInnerSymmetry      =   muonBarrelExtension->inner_symmetry;
    m_caloHitCreatorSettings.m_hCalEndCapOuterR             =   hCalEndcapExtension->extent[1]/dd4hep::mm;
    m_caloHitCreatorSettings.m_hCalEndCapOuterZ             =   hCalEndcapExtension->extent[3]/dd4hep::mm;
    m_caloHitCreatorSettings.m_hCalBarrelOuterR             =   hCalBarrelExtension->extent[1]/dd4hep::mm;
    m_caloHitCreatorSettings.m_hCalBarrelOuterPhi0          =   hCalBarrelExtension->outer_phi0/dd4hep::rad;
    m_caloHitCreatorSettings.m_hCalBarrelOuterSymmetry      =   hCalBarrelExtension->outer_symmetry;
    m_caloHitCreatorSettings.m_hCalEndCapInnerSymmetryOrder =   hCalEndcapExtension->inner_symmetry;;
    m_caloHitCreatorSettings.m_hCalEndCapInnerPhiCoordinate =   hCalEndcapExtension->inner_phi0/dd4hep::rad;;
    
    // Get the magnetic field
    DD4hep::Geometry::LCDD& lcdd = DD4hep::Geometry::LCDD::getInstance();
    const double position[3]={0,0,0}; // position to calculate magnetic field at (the origin in this case)
    double magneticFieldVector[3]={0,0,0}; // initialise object to hold magnetic field
    lcdd.field().magneticField(position,magneticFieldVector); // get the magnetic field vector from DD4hep
    
    m_settings.m_innerBField = magneticFieldVector[2]/dd4hep::tesla; // z component at (0,0,0)
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
