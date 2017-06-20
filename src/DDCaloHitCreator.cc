/**
 *  @file   DDMarlinPandora/src/DDCaloHitCreator.cc
 * 
 *  @brief  Implementation of the calo hit creator class.
 * 
 *  $Log: $
 */

#include "DDCaloHitCreator.h"

#include "marlin/Global.h"
#include "marlin/Processor.h"

#include "UTIL/CellIDDecoder.h"


#include <DD4hep/DD4hepUnits.h>
#include <DD4hep/DetType.h>
#include <DD4hep/DetectorSelector.h>
#include <DD4hep/Detector.h>

#include <algorithm>
#include <cmath>
#include <limits>

//forward declarations. See in DDPandoraPFANewProcessor.cc

// dd4hep::rec::LayeredCalorimeterData * getExtension(std::string detectorName);
dd4hep::rec::LayeredCalorimeterData * getExtension(unsigned int includeFlag, unsigned int excludeFlag=0);

// double getCoilOuterR();

///FIXME: HANDLE PROBLEM WHEN EXTENSION IS MISSING
DDCaloHitCreator::DDCaloHitCreator(const Settings &settings, const pandora::Pandora *const pPandora) :
    m_settings(settings),
    m_pandora(*pPandora),
    m_hCalBarrelLayerThickness(0.f),
    m_hCalEndCapLayerThickness(0.f),
    m_calorimeterHitVector(0),
    m_volumeManager()
{
    
    const std::vector<dd4hep::rec::LayeredCalorimeterStruct::Layer>& barrelLayers= getExtension(( dd4hep::DetType::CALORIMETER | dd4hep::DetType::HADRONIC | dd4hep::DetType::BARREL), ( dd4hep::DetType::AUXILIARY  |  dd4hep::DetType::FORWARD ))->layers;
    
    const std::vector<dd4hep::rec::LayeredCalorimeterStruct::Layer>& endcapLayers= getExtension(( dd4hep::DetType::CALORIMETER | dd4hep::DetType::HADRONIC | dd4hep::DetType::ENDCAP), ( dd4hep::DetType::AUXILIARY  |  dd4hep::DetType::FORWARD ))->layers;
    
    ///Take thicknesses from last layer (was like that before with gear)
    m_hCalEndCapLayerThickness =(endcapLayers.back().inner_thickness+endcapLayers.back().outer_thickness)/dd4hep::mm;
    m_hCalBarrelLayerThickness =(barrelLayers.back().inner_thickness+barrelLayers.back().outer_thickness)/dd4hep::mm;

    if ((m_hCalEndCapLayerThickness < std::numeric_limits<float>::epsilon()) || (m_hCalBarrelLayerThickness < std::numeric_limits<float>::epsilon()))
        throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);

   dd4hep::Detector& theDetector = dd4hep::Detector::getInstance();
   m_volumeManager = theDetector.volumeManager();
   if( not m_volumeManager.isValid() ){
     theDetector.apply("DD4hepVolumeManager",0,0);
     m_volumeManager = theDetector.volumeManager();
   }

}

//------------------------------------------------------------------------------------------------------------------------------------------

DDCaloHitCreator::~DDCaloHitCreator()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode DDCaloHitCreator::CreateCaloHits(const EVENT::LCEvent *const pLCEvent)
{

  //fg: there cannot be any reasonable default for this string - so we set it to sth. that will cause an exception in case 
  //    the cellID encoding string is not in the collection: 
  UTIL::CellIDDecoder<CalorimeterHit>::setDefaultEncoding("undefined_cellID_encoding:100");
  
    PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, this->CreateECalCaloHits(pLCEvent));
    PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, this->CreateHCalCaloHits(pLCEvent));
    PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, this->CreateMuonCaloHits(pLCEvent));
    PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, this->CreateLCalCaloHits(pLCEvent));
    PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, this->CreateLHCalCaloHits(pLCEvent));

    return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode DDCaloHitCreator::CreateECalCaloHits(const EVENT::LCEvent *const pLCEvent)
{
    for (StringVector::const_iterator iter = m_settings.m_eCalCaloHitCollections.begin(), iterEnd = m_settings.m_eCalCaloHitCollections.end();
        iter != iterEnd; ++iter)
    {
        try
        {
            const EVENT::LCCollection *pCaloHitCollection = pLCEvent->getCollection(*iter);
            const int nElements(pCaloHitCollection->getNumberOfElements());

            if (0 == nElements)
                continue;

	    streamlog_out( DEBUG1 ) << "Creating " << *iter << " hits" << std::endl;

            const std::vector<dd4hep::rec::LayeredCalorimeterStruct::Layer>& barrelLayers= getExtension( ( dd4hep::DetType::CALORIMETER | dd4hep::DetType::ELECTROMAGNETIC | dd4hep::DetType::BARREL), ( dd4hep::DetType::AUXILIARY  |  dd4hep::DetType::FORWARD ) )->layers;
            const std::vector<dd4hep::rec::LayeredCalorimeterStruct::Layer>& endcapLayers= getExtension( ( dd4hep::DetType::CALORIMETER | dd4hep::DetType::ELECTROMAGNETIC | dd4hep::DetType::ENDCAP), ( dd4hep::DetType::AUXILIARY  |  dd4hep::DetType::FORWARD ) )->layers;
            

            UTIL::CellIDDecoder<CalorimeterHit> cellIdDecoder(pCaloHitCollection);
            const std::string layerCodingString(pCaloHitCollection->getParameters().getStringVal(LCIO::CellIDEncoding));
            std::string layerCoding("layer");

            for (int i = 0; i < nElements; ++i)
            {
                try
                {
                    EVENT::CalorimeterHit *pCaloHit = dynamic_cast<CalorimeterHit*>(pCaloHitCollection->getElementAt(i));

                    if (NULL == pCaloHit)
                        throw EVENT::Exception("Collection type mismatch");

                    float eCalToMip(m_settings.m_eCalToMip), eCalMipThreshold(m_settings.m_eCalMipThreshold), eCalToEMGeV(m_settings.m_eCalToEMGeV),
                        eCalToHadGeVBarrel(m_settings.m_eCalToHadGeVBarrel), eCalToHadGeVEndCap(m_settings.m_eCalToHadGeVEndCap);

                    // Hybrid ECAL including pure ScECAL.
                    if (m_settings.m_useEcalScLayers)
                    {
                        std::string collectionName(*iter);
                        std::transform(collectionName.begin(), collectionName.end(), collectionName.begin(), ::tolower);
                        layerCoding = "layer";

                        if (collectionName.find("ecal", 0) == std::string::npos)
                            streamlog_out(MESSAGE) << "WARNING: mismatching hybrid Ecal collection name. " << collectionName << std::endl;

                        if (collectionName.find("si", 0) != std::string::npos)
                        {
                             eCalToMip = m_settings.m_eCalSiToMip;
                             eCalMipThreshold = m_settings.m_eCalSiMipThreshold;
                             eCalToEMGeV = m_settings.m_eCalSiToEMGeV;
                             eCalToHadGeVBarrel = m_settings.m_eCalSiToHadGeVBarrel;
                             eCalToHadGeVEndCap = m_settings.m_eCalSiToHadGeVEndCap;
                        }
                        else if (collectionName.find("sc", 0) != std::string::npos)
                        {
                             eCalToMip = m_settings.m_eCalScToMip;
                             eCalMipThreshold = m_settings.m_eCalScMipThreshold;
                             eCalToEMGeV = m_settings.m_eCalScToEMGeV;
                             eCalToHadGeVBarrel = m_settings.m_eCalScToHadGeVBarrel;
                             eCalToHadGeVEndCap = m_settings.m_eCalScToHadGeVEndCap;
                        }
                    }

                    PandoraApi::CaloHit::Parameters caloHitParameters;
                    caloHitParameters.m_hitType = pandora::ECAL;
                    caloHitParameters.m_isDigital = false;
                    caloHitParameters.m_layer = cellIdDecoder(pCaloHit)[layerCoding.c_str()];
                    caloHitParameters.m_isInOuterSamplingLayer = false;
                    this->GetCommonCaloHitProperties(pCaloHit, caloHitParameters);

                    float absorberCorrection(1.);

                    //FIXME: Why is this used to get barrel or endcap??? use SystemID if available
                    if (std::fabs(pCaloHit->getPosition()[2]) < m_settings.m_eCalBarrelOuterZ)
                    {

                      streamlog_out( DEBUG6 ) << "IDS " << *iter
					      << std::setw(15) << pCaloHit->getCellID0()
					      << std::setw(15) << pCaloHit->getPosition()[0]
					      << std::setw(15) << pCaloHit->getPosition()[1]
					      << std::setw(15) << pCaloHit->getPosition()[2]
					      << std::setw(5) << cellIdDecoder(pCaloHit)["system"]
					      << std::setw(5) << cellIdDecoder(pCaloHit)["module"]
					      << std::setw(5) << cellIdDecoder(pCaloHit)["stave"]
					      << std::setw(5) << cellIdDecoder(pCaloHit)["layer"]
					      << std::setw(5) << cellIdDecoder(pCaloHit)["x"]
					      << std::setw(5) << cellIdDecoder(pCaloHit)["y"]
					      << std::endl;

                      this->GetBarrelCaloHitProperties(pCaloHit, barrelLayers, m_settings.m_eCalBarrelInnerSymmetry, caloHitParameters, m_settings.m_eCalBarrelNormalVector, absorberCorrection);

                        caloHitParameters.m_hadronicEnergy = eCalToHadGeVBarrel * pCaloHit->getEnergy();
                    }
                    else
                    {
                        this->GetEndCapCaloHitProperties(pCaloHit, endcapLayers, caloHitParameters, absorberCorrection);
                        caloHitParameters.m_hadronicEnergy = eCalToHadGeVEndCap * pCaloHit->getEnergy();
                    }

                    caloHitParameters.m_mipEquivalentEnergy = pCaloHit->getEnergy() * eCalToMip * absorberCorrection;

                    if (caloHitParameters.m_mipEquivalentEnergy.Get() < eCalMipThreshold)
                        continue;

                    caloHitParameters.m_electromagneticEnergy = eCalToEMGeV * pCaloHit->getEnergy();

                    // ATTN If using strip splitting, must correct cell sizes for use in PFA to minimum of strip width and strip length
                    if (m_settings.m_stripSplittingOn)
                    {
                        const float splitCellSize(std::min(caloHitParameters.m_cellSize0.Get(), caloHitParameters.m_cellSize1.Get()));
                        caloHitParameters.m_cellSize0 = splitCellSize;
                        caloHitParameters.m_cellSize1 = splitCellSize;
                    }

                    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::CaloHit::Create(m_pandora, caloHitParameters));
                    m_calorimeterHitVector.push_back(pCaloHit);

                }
                catch (pandora::StatusCodeException &statusCodeException)
                {
                    streamlog_out(ERROR) << "Failed to extract ecal calo hit from "
					 << *iter <<  ": "
					 << statusCodeException.ToString() << std::endl;
                }
                catch (EVENT::Exception &exception)
                {
                    streamlog_out(WARNING) << "Failed to extract ecal calo hit from "
					   << *iter <<  ": "
					   << exception.what() << std::endl;
                }
            }
        }
        catch (EVENT::Exception &exception)
        {
            streamlog_out(MESSAGE) << "Failed to extract ecal calo hit collection: " << *iter << ", " << exception.what() << std::endl;
        }
    }

    return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode DDCaloHitCreator::CreateHCalCaloHits(const EVENT::LCEvent *const pLCEvent)
{
    for (StringVector::const_iterator iter = m_settings.m_hCalCaloHitCollections.begin(), iterEnd = m_settings.m_hCalCaloHitCollections.end();
        iter != iterEnd; ++iter)
    {
        try
        {
            const EVENT::LCCollection *pCaloHitCollection = pLCEvent->getCollection(*iter);
            const int nElements(pCaloHitCollection->getNumberOfElements());
            UTIL::CellIDDecoder<CalorimeterHit> cellIdDecoder(pCaloHitCollection);
            
            if (0 == nElements)
                continue;
            
	    streamlog_out( DEBUG1 ) << "Creating " << *iter << " hits" << std::endl;

            const std::vector<dd4hep::rec::LayeredCalorimeterStruct::Layer>& barrelLayers= getExtension( ( dd4hep::DetType::CALORIMETER | dd4hep::DetType::HADRONIC | dd4hep::DetType::BARREL), ( dd4hep::DetType::AUXILIARY  |  dd4hep::DetType::FORWARD ) )->layers;
            const std::vector<dd4hep::rec::LayeredCalorimeterStruct::Layer>& endcapLayers= getExtension( ( dd4hep::DetType::CALORIMETER | dd4hep::DetType::HADRONIC| dd4hep::DetType::ENDCAP), ( dd4hep::DetType::AUXILIARY )  |  dd4hep::DetType::FORWARD )->layers;
            
            
            
            const std::string layerCoding("layer");

            for (int i = 0; i < nElements; ++i)
            {
                try
                {
                    EVENT::CalorimeterHit *pCaloHit = dynamic_cast<CalorimeterHit*>(pCaloHitCollection->getElementAt(i));

                    if (NULL == pCaloHit)
                        throw EVENT::Exception("Collection type mismatch");

                    PandoraApi::CaloHit::Parameters caloHitParameters;
                    caloHitParameters.m_hitType = pandora::HCAL;
                    caloHitParameters.m_isDigital = false;
                    caloHitParameters.m_layer = cellIdDecoder(pCaloHit)[layerCoding.c_str()];
                    caloHitParameters.m_isInOuterSamplingLayer = (this->GetNLayersFromEdge(pCaloHit) <= m_settings.m_nOuterSamplingLayers);
                    this->GetCommonCaloHitProperties(pCaloHit, caloHitParameters);

                    float absorberCorrection(1.);

                    if (std::fabs(pCaloHit->getPosition()[2]) < m_settings.m_hCalBarrelOuterZ)
                    {

                      this->GetBarrelCaloHitProperties(pCaloHit, barrelLayers, m_settings.m_hCalBarrelInnerSymmetry, caloHitParameters, m_settings.m_hCalBarrelNormalVector, absorberCorrection);
                    }
                    else
                    {
                        this->GetEndCapCaloHitProperties(pCaloHit, endcapLayers, caloHitParameters, absorberCorrection);
                    }

                    caloHitParameters.m_mipEquivalentEnergy = pCaloHit->getEnergy() * m_settings.m_hCalToMip * absorberCorrection;

                    if (caloHitParameters.m_mipEquivalentEnergy.Get() < m_settings.m_hCalMipThreshold)
                        continue;

                    caloHitParameters.m_hadronicEnergy = std::min(m_settings.m_hCalToHadGeV * pCaloHit->getEnergy(), m_settings.m_maxHCalHitHadronicEnergy);
                    caloHitParameters.m_electromagneticEnergy = m_settings.m_hCalToEMGeV * pCaloHit->getEnergy();

                    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::CaloHit::Create(m_pandora, caloHitParameters));
                    m_calorimeterHitVector.push_back(pCaloHit);
                }
                catch (pandora::StatusCodeException &statusCodeException)
                {
                    streamlog_out(ERROR) << "Failed to extract hcal calo hit: " << statusCodeException.ToString() << std::endl;
                }
                catch (EVENT::Exception &exception)
                {
                    streamlog_out(WARNING) << "Failed to extract hcal calo hit: " << exception.what() << std::endl;
                }
            }
        }
        catch (EVENT::Exception &exception)
        {
            streamlog_out(MESSAGE) << "Failed to extract hcal calo hit collection: " << *iter << ", " << exception.what() << std::endl;
        }
    }

    return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode DDCaloHitCreator::CreateMuonCaloHits(const EVENT::LCEvent *const pLCEvent)
{
    for (StringVector::const_iterator iter = m_settings.m_muonCaloHitCollections.begin(), iterEnd = m_settings.m_muonCaloHitCollections.end();
        iter != iterEnd; ++iter)
    {
        try
        {
            const EVENT::LCCollection *pCaloHitCollection = pLCEvent->getCollection(*iter);
            const int nElements(pCaloHitCollection->getNumberOfElements());

            if (0 == nElements)
                continue;

	    streamlog_out( DEBUG1 ) << "Creating " << *iter << " hits" << std::endl;

            const std::vector<dd4hep::rec::LayeredCalorimeterStruct::Layer>& barrelLayers= getExtension(( dd4hep::DetType::CALORIMETER | dd4hep::DetType::MUON| dd4hep::DetType::BARREL), ( dd4hep::DetType::AUXILIARY  |  dd4hep::DetType::FORWARD  ))->layers;
            const std::vector<dd4hep::rec::LayeredCalorimeterStruct::Layer>& endcapLayers= getExtension(( dd4hep::DetType::CALORIMETER | dd4hep::DetType::MUON| dd4hep::DetType::ENDCAP), ( dd4hep::DetType::AUXILIARY  |  dd4hep::DetType::FORWARD  ))->layers;
            ///FIXME: WHAT ABOUT MORE MUON SYSTEMS?
            // const std::vector<dd4hep::rec::LayeredCalorimeterStruct::Layer>& plugLayers= getExtension(( dd4hep::DetType::CALORIMETER | dd4hep::DetType::HADRONIC | dd4hep::DetType::AUXILIARY ))->layers;
            UTIL::CellIDDecoder<CalorimeterHit> cellIdDecoder(pCaloHitCollection);
            const std::string layerCodingString(pCaloHitCollection->getParameters().getStringVal(LCIO::CellIDEncoding));
            const std::string layerCoding("layer");

            for (int i = 0; i < nElements; ++i)
            {
                try
                {
                    EVENT::CalorimeterHit *pCaloHit = dynamic_cast<CalorimeterHit*>(pCaloHitCollection->getElementAt(i));

                    if (NULL == pCaloHit)
                        throw EVENT::Exception("Collection type mismatch");

                    PandoraApi::CaloHit::Parameters caloHitParameters;
                    caloHitParameters.m_hitType = pandora::MUON;
                    caloHitParameters.m_layer = cellIdDecoder(pCaloHit)[layerCoding.c_str()];
                    caloHitParameters.m_isInOuterSamplingLayer = true;
                    this->GetCommonCaloHitProperties(pCaloHit, caloHitParameters);

                    const float radius(std::sqrt(pCaloHit->getPosition()[0] * pCaloHit->getPosition()[0] +
                        pCaloHit->getPosition()[1] * pCaloHit->getPosition()[1]));

                    const bool isWithinCoil(radius < m_settings.m_coilOuterR);
                    const bool isInBarrelRegion(std::fabs(pCaloHit->getPosition()[2]) < m_settings.m_muonBarrelOuterZ);

                    float absorberCorrection(1.);

                    if (isInBarrelRegion && isWithinCoil)
                    {
                      std::cout<<"BIG WARNING: CANNOT HANDLE PLUG HITS (no plug in CLIC model), DO NOTHING!"<<std::endl;
//                         this->GetEndCapCaloHitProperties(pCaloHit, plugLayers, caloHitParameters, absorberCorrection);
                    }
                    else if (isInBarrelRegion)
                    {

                      this->GetBarrelCaloHitProperties(pCaloHit, barrelLayers, m_settings.m_muonBarrelInnerSymmetry, caloHitParameters, m_settings.m_muonBarrelNormalVector, absorberCorrection);
                    }
                    else
                    {
                        this->GetEndCapCaloHitProperties(pCaloHit, endcapLayers, caloHitParameters, absorberCorrection);
                    }

                    if (m_settings.m_muonDigitalHits > 0)
                    {
                        caloHitParameters.m_isDigital = true;
                        caloHitParameters.m_inputEnergy = m_settings.m_muonHitEnergy;
                        caloHitParameters.m_hadronicEnergy = m_settings.m_muonHitEnergy;
                        caloHitParameters.m_electromagneticEnergy = m_settings.m_muonHitEnergy;
                        caloHitParameters.m_mipEquivalentEnergy = 1.f;
                    }
                    else
                    {
                        caloHitParameters.m_isDigital = false;
                        caloHitParameters.m_inputEnergy = pCaloHit->getEnergy();
                        caloHitParameters.m_hadronicEnergy = pCaloHit->getEnergy();
                        caloHitParameters.m_electromagneticEnergy = pCaloHit->getEnergy();
                        caloHitParameters.m_mipEquivalentEnergy = pCaloHit->getEnergy() * m_settings.m_muonToMip;
                    }

                    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::CaloHit::Create(m_pandora, caloHitParameters));
                    m_calorimeterHitVector.push_back(pCaloHit);
                }
                catch (pandora::StatusCodeException &statusCodeException)
                {
                    streamlog_out(ERROR) << "Failed to extract muon hit: " << statusCodeException.ToString() << std::endl;
                }
                catch (EVENT::Exception &exception)
                {
                    streamlog_out(WARNING) << "Failed to extract muon hit: " << exception.what() << std::endl;
                }
            }
        }
        catch (EVENT::Exception &exception)
        {
            streamlog_out(MESSAGE) << "Failed to extract muon hit collection: " << *iter << ", " << exception.what() << std::endl;
        }
    }

    return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode DDCaloHitCreator::CreateLCalCaloHits(const EVENT::LCEvent *const pLCEvent)
{
    for (StringVector::const_iterator iter = m_settings.m_lCalCaloHitCollections.begin(), iterEnd = m_settings.m_lCalCaloHitCollections.end();
        iter != iterEnd; ++iter)
    {
        try
        {
            const EVENT::LCCollection *pCaloHitCollection = pLCEvent->getCollection(*iter);
            const int nElements(pCaloHitCollection->getNumberOfElements());

            if (0 == nElements)
                continue;

	    streamlog_out( DEBUG1 ) << "Creating " << *iter << " hits" << std::endl;

            ///FIXME: WHAT ABOUT OTHER ECALS?
            const std::vector<dd4hep::rec::LayeredCalorimeterStruct::Layer>& endcapLayers= getExtension( dd4hep::DetType::CALORIMETER |  dd4hep::DetType::ENDCAP  | dd4hep::DetType::ELECTROMAGNETIC |  dd4hep::DetType::FORWARD ,  dd4hep::DetType::AUXILIARY )->layers;
        

            UTIL::CellIDDecoder<CalorimeterHit> cellIdDecoder(pCaloHitCollection);
            const std::string layerCodingString(pCaloHitCollection->getParameters().getStringVal(LCIO::CellIDEncoding));
            const std::string layerCoding("layer");

            for (int i = 0; i < nElements; ++i)
            {
                try
                {
                    EVENT::CalorimeterHit *pCaloHit = dynamic_cast<CalorimeterHit*>(pCaloHitCollection->getElementAt(i));

                    if (NULL == pCaloHit)
                        throw EVENT::Exception("Collection type mismatch");

                    PandoraApi::CaloHit::Parameters caloHitParameters;
                    caloHitParameters.m_hitType = pandora::ECAL;
                    caloHitParameters.m_isDigital = false;
                    caloHitParameters.m_layer = cellIdDecoder(pCaloHit)[layerCoding.c_str()];
                    caloHitParameters.m_isInOuterSamplingLayer = false;
                    this->GetCommonCaloHitProperties(pCaloHit, caloHitParameters);

                    float absorberCorrection(1.);
                    this->GetEndCapCaloHitProperties(pCaloHit, endcapLayers, caloHitParameters, absorberCorrection);

                    caloHitParameters.m_mipEquivalentEnergy = pCaloHit->getEnergy() * m_settings.m_eCalToMip * absorberCorrection;

                    if (caloHitParameters.m_mipEquivalentEnergy.Get() < m_settings.m_eCalMipThreshold)
                        continue;

                    caloHitParameters.m_electromagneticEnergy = m_settings.m_eCalToEMGeV * pCaloHit->getEnergy();
                    caloHitParameters.m_hadronicEnergy = m_settings.m_eCalToHadGeVEndCap * pCaloHit->getEnergy();

                    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::CaloHit::Create(m_pandora, caloHitParameters));
                    m_calorimeterHitVector.push_back(pCaloHit);
                }
                catch (pandora::StatusCodeException &statusCodeException)
                {
                    streamlog_out(ERROR) << "Failed to extract lcal calo hit: " << statusCodeException.ToString() << std::endl;
                }
                catch (EVENT::Exception &exception)
                {
                    streamlog_out(WARNING) << "Failed to extract lcal calo hit: " << exception.what() << std::endl;
                }
            }
        }
        catch (EVENT::Exception &exception)
        {
            streamlog_out(MESSAGE) << "Failed to extract lcal calo hit collection: " << *iter << ", " << exception.what() << std::endl;
        }
    }

    return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode DDCaloHitCreator::CreateLHCalCaloHits(const EVENT::LCEvent *const pLCEvent)
{
    for (StringVector::const_iterator iter = m_settings.m_lHCalCaloHitCollections.begin(), iterEnd = m_settings.m_lHCalCaloHitCollections.end();
        iter != iterEnd; ++iter)
    {
        try
        {
            const EVENT::LCCollection *pCaloHitCollection = pLCEvent->getCollection(*iter);
            const int nElements(pCaloHitCollection->getNumberOfElements());

            if (0 == nElements)
                continue;

	    streamlog_out( DEBUG1 ) << "Creating " << *iter << " hits" << std::endl;

            ///FIXME! WHAT ABOUT MORE HCALS?
            const std::vector<dd4hep::rec::LayeredCalorimeterStruct::Layer>& endcapLayers= getExtension(dd4hep::DetType::CALORIMETER |  dd4hep::DetType::ENDCAP  | dd4hep::DetType::HADRONIC|  dd4hep::DetType::FORWARD)->layers;
            
            UTIL::CellIDDecoder<CalorimeterHit> cellIdDecoder(pCaloHitCollection);
            const std::string layerCodingString(pCaloHitCollection->getParameters().getStringVal(LCIO::CellIDEncoding));
            const std::string layerCoding("layer");

            for (int i = 0; i < nElements; ++i)
            {
                try
                {
                    EVENT::CalorimeterHit *pCaloHit = dynamic_cast<CalorimeterHit*>(pCaloHitCollection->getElementAt(i));

                    if (NULL == pCaloHit)
                        throw EVENT::Exception("Collection type mismatch");

                    PandoraApi::CaloHit::Parameters caloHitParameters;
                    caloHitParameters.m_hitType = pandora::HCAL;
                    caloHitParameters.m_isDigital = false;
                    caloHitParameters.m_layer = cellIdDecoder(pCaloHit)[layerCoding.c_str()];
                    caloHitParameters.m_isInOuterSamplingLayer = (this->GetNLayersFromEdge(pCaloHit) <= m_settings.m_nOuterSamplingLayers);
                    this->GetCommonCaloHitProperties(pCaloHit, caloHitParameters);

                    float absorberCorrection(1.);
                    this->GetEndCapCaloHitProperties(pCaloHit, endcapLayers, caloHitParameters, absorberCorrection);

                    caloHitParameters.m_mipEquivalentEnergy = pCaloHit->getEnergy() * m_settings.m_hCalToMip * absorberCorrection;

                    if (caloHitParameters.m_mipEquivalentEnergy.Get() < m_settings.m_hCalMipThreshold)
                        continue;

                    caloHitParameters.m_hadronicEnergy = std::min(m_settings.m_hCalToHadGeV * pCaloHit->getEnergy(), m_settings.m_maxHCalHitHadronicEnergy);
                    caloHitParameters.m_electromagneticEnergy = m_settings.m_hCalToEMGeV * pCaloHit->getEnergy();

                    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::CaloHit::Create(m_pandora, caloHitParameters));
                    m_calorimeterHitVector.push_back(pCaloHit);
                }
                catch (pandora::StatusCodeException &statusCodeException)
                {
                    streamlog_out(ERROR) << "Failed to extract lhcal calo hit: " << statusCodeException.ToString() << std::endl;
                }
                catch (EVENT::Exception &exception)
                {
                    streamlog_out(WARNING) << "Failed to extract lhcal calo hit: " << exception.what() << std::endl;
                }
            }
        }
        catch (EVENT::Exception &exception)
        {
            streamlog_out(MESSAGE) << "Failed to extract lhcal calo hit collection: " << *iter << ", " << exception.what() << std::endl;
        }
    }

    return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DDCaloHitCreator::GetCommonCaloHitProperties(const EVENT::CalorimeterHit *const pCaloHit, PandoraApi::CaloHit::Parameters &caloHitParameters) const
{
    const float *pCaloHitPosition(pCaloHit->getPosition());
    const pandora::CartesianVector positionVector(pCaloHitPosition[0], pCaloHitPosition[1], pCaloHitPosition[2]);

    caloHitParameters.m_cellGeometry = pandora::RECTANGULAR;
    caloHitParameters.m_positionVector = positionVector;
    caloHitParameters.m_expectedDirection = positionVector.GetUnitVector();
    caloHitParameters.m_pParentAddress = (void*)pCaloHit;
    caloHitParameters.m_inputEnergy = pCaloHit->getEnergy();
    caloHitParameters.m_time = pCaloHit->getTime();
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DDCaloHitCreator::GetEndCapCaloHitProperties(const EVENT::CalorimeterHit *const pCaloHit, const std::vector<dd4hep::rec::LayeredCalorimeterStruct::Layer> &layers,
    PandoraApi::CaloHit::Parameters &caloHitParameters, float &absorberCorrection) const
{
    caloHitParameters.m_hitRegion = pandora::ENDCAP;

    //FIXME! WHAT DO WE DO HERE?
    const int physicalLayer(std::min(static_cast<int>(caloHitParameters.m_layer.Get()), static_cast<int>(layers.size()-1)));
    caloHitParameters.m_cellSize0 = layers[physicalLayer].cellSize0/dd4hep::mm;
    caloHitParameters.m_cellSize1 = layers[physicalLayer].cellSize1/dd4hep::mm;
    
    double thickness = (layers[physicalLayer].inner_thickness+layers[physicalLayer].sensitive_thickness/2.0)/dd4hep::mm;
    double nRadLengths = layers[physicalLayer].inner_nRadiationLengths;
    double nIntLengths = layers[physicalLayer].inner_nInteractionLengths;
    double layerAbsorberThickness = (layers[physicalLayer].inner_thickness-layers[physicalLayer].sensitive_thickness/2.0)/dd4hep::mm;

    if(physicalLayer>0){
        thickness += (layers[physicalLayer-1].outer_thickness -layers[physicalLayer].sensitive_thickness/2.0)/dd4hep::mm;
        nRadLengths += layers[physicalLayer-1].outer_nRadiationLengths;
        nIntLengths += layers[physicalLayer-1].outer_nInteractionLengths;
        layerAbsorberThickness += (layers[physicalLayer-1].outer_thickness -layers[physicalLayer].sensitive_thickness/2.0)/dd4hep::mm;

    }
    
    caloHitParameters.m_cellThickness = thickness;
    caloHitParameters.m_nCellRadiationLengths = nRadLengths;
    caloHitParameters.m_nCellInteractionLengths = nIntLengths;
    
    if (caloHitParameters.m_nCellRadiationLengths.Get() < std::numeric_limits<float>::epsilon() || caloHitParameters.m_nCellInteractionLengths.Get() < std::numeric_limits<float>::epsilon())
    {
        streamlog_out(WARNING) << "CaloHitCreator::GetEndCapCaloHitProperties Calo hit has 0 radiation length or interaction length: \
            not creating a Pandora calo hit." << std::endl;
        throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);
    }

    
    //FIXME! do we need this?
    absorberCorrection = 1.;
    for (unsigned int i = 0, iMax = layers.size(); i < iMax; ++i)
    {
        float absorberThickness((layers[i].inner_thickness - layers[i].sensitive_thickness/2.0 )/dd4hep::mm);
        
        if (i>0)
            absorberThickness += (layers[i-1].outer_thickness - layers[i-1].sensitive_thickness/2.0)/dd4hep::mm;

        if (absorberThickness < std::numeric_limits<float>::epsilon())
            continue;

        if (layerAbsorberThickness > std::numeric_limits<float>::epsilon())
            absorberCorrection = absorberThickness / layerAbsorberThickness;

        break;
    }

    caloHitParameters.m_cellNormalVector = (pCaloHit->getPosition()[2] > 0) ? pandora::CartesianVector(0, 0, 1) :
        pandora::CartesianVector(0, 0, -1);
        
//     streamlog_out(DEBUG) <<" GetEndCapCaloHitProperties: physLayer: "<<physicalLayer <<" layer: "<<caloHitParameters.m_layer.Get()<<" nX0: "<<    caloHitParameters.m_nCellRadiationLengths.Get() <<" nLambdaI: "<<    caloHitParameters.m_nCellInteractionLengths.Get()<<" thickness: "<<caloHitParameters.m_cellThickness.Get()<<std::endl;
    
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DDCaloHitCreator::GetBarrelCaloHitProperties( const EVENT::CalorimeterHit *const pCaloHit,
                                                   const std::vector<dd4hep::rec::LayeredCalorimeterStruct::Layer> &layers,
                                                   unsigned int barrelSymmetryOrder,
                                                   PandoraApi::CaloHit::Parameters &caloHitParameters,
                                                   FloatVector const& normalVector,
                                                   float &absorberCorrection ) const
{
    caloHitParameters.m_hitRegion = pandora::BARREL;

    //FIXME! WHAT DO WE DO HERE?
    const int physicalLayer(std::min(static_cast<int>(caloHitParameters.m_layer.Get()), static_cast<int>(layers.size()-1)));
    caloHitParameters.m_cellSize0 = layers[physicalLayer].cellSize0/dd4hep::mm;
    caloHitParameters.m_cellSize1 = layers[physicalLayer].cellSize1/dd4hep::mm;
    
    double thickness = (layers[physicalLayer].inner_thickness+layers[physicalLayer].sensitive_thickness/2.0)/dd4hep::mm;
    double nRadLengths = layers[physicalLayer].inner_nRadiationLengths;
    double nIntLengths = layers[physicalLayer].inner_nInteractionLengths;

    double layerAbsorberThickness = (layers[physicalLayer].inner_thickness-layers[physicalLayer].sensitive_thickness/2.0)/dd4hep::mm;
    if(physicalLayer>0){
        thickness += (layers[physicalLayer-1].outer_thickness -layers[physicalLayer].sensitive_thickness/2.0)/dd4hep::mm;
        nRadLengths += layers[physicalLayer-1].outer_nRadiationLengths;
        nIntLengths += layers[physicalLayer-1].outer_nInteractionLengths;
        layerAbsorberThickness += (layers[physicalLayer-1].outer_thickness -layers[physicalLayer].sensitive_thickness/2.0)/dd4hep::mm;
    }
    
    caloHitParameters.m_cellThickness = thickness;
    caloHitParameters.m_nCellRadiationLengths = nRadLengths;
    caloHitParameters.m_nCellInteractionLengths = nIntLengths;

    if (caloHitParameters.m_nCellRadiationLengths.Get() < std::numeric_limits<float>::epsilon() || caloHitParameters.m_nCellInteractionLengths.Get() < std::numeric_limits<float>::epsilon())
    {
        streamlog_out(WARNING) << "CaloHitCreator::GetBarrelCaloHitProperties Calo hit has 0 radiation length or interaction length: \
            not creating a Pandora calo hit." << std::endl;
        throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);
    }

    //FIXME! do we need this?
    absorberCorrection = 1.;
    for (unsigned int i = 0, iMax = layers.size(); i < iMax; ++i)
    {
        float absorberThickness((layers[i].inner_thickness - layers[i].sensitive_thickness/2.0 )/dd4hep::mm);
        
        if (i>0)
            absorberThickness += (layers[i-1].outer_thickness - layers[i-1].sensitive_thickness/2.0)/dd4hep::mm;

        if (absorberThickness < std::numeric_limits<float>::epsilon())
            continue;

        if (layerAbsorberThickness > std::numeric_limits<float>::epsilon())
            absorberCorrection = absorberThickness / layerAbsorberThickness;

        break;
    }

    if (barrelSymmetryOrder > 2)
    {

      if ( pCaloHit->getCellID0() != 0 ) {

        auto staveDetElement = m_volumeManager.lookupDetElement( pCaloHit->getCellID0() );
        dd4hep::Position local1(0.0, 0.0, 0.0);
        dd4hep::Position local2(normalVector[0],normalVector[1],normalVector[2]);
        dd4hep::Position global1(0.0, 0.0, 0.0);
        dd4hep::Position global2(0.0, 0.0, 0.0);
        staveDetElement.nominal().localToWorld( local1, global1 );
        staveDetElement.nominal().localToWorld( local2, global2 );
        dd4hep::Position normal( global2-global1 );

        streamlog_out(DEBUG6) << "   detelement: " << staveDetElement.name()
			      << "   parent: " << staveDetElement.parent().name()
			      << "   grandparent: " << staveDetElement.parent().parent().name()
			      << "   cellID: " << pCaloHit->getCellID0()
                              << "   PhiLoc:"  << atan2( global1.y(), global1.x() )*180/M_PI
                              << "   PhiNor:"  << atan2( normal.y(), normal.x() )*180/M_PI
                              << " normal vector "
                              << std::setw(15) << normal.x()
                              << std::setw(15) << normal.y()
                              << std::setw(15) << normal.z()
                              << std::endl;

        caloHitParameters.m_cellNormalVector = pandora::CartesianVector( normal.x(), normal.y(), normal.z() );
      } else {
	const double phi = atan2( pCaloHit->getPosition()[1], pCaloHit->getPosition()[0] );

	streamlog_out( WARNING )  << "This hit does not have any cellIDs set, will use phi-direction for normal vector "
				  << " phi:" << std::setw(15) << phi*180/M_PI
				  << std::endl;

        caloHitParameters.m_cellNormalVector = pandora::CartesianVector( std::cos(phi), std::sin(phi) , 0.0 );
      }



    }
    else
    {
        const float *pCaloHitPosition( pCaloHit->getPosition() );
        const float phi = std::atan2( pCaloHitPosition[1], pCaloHitPosition[0] );
        caloHitParameters.m_cellNormalVector = pandora::CartesianVector(std::cos(phi), std::sin(phi), 0);
    }
    
//     streamlog_out(DEBUG)<<" GetBarrelCaloHitProperties: physLayer: "<<physicalLayer <<" layer: "<<caloHitParameters.m_layer.Get()<<" nX0: "<<    caloHitParameters.m_nCellRadiationLengths.Get() <<" nLambdaI: "<<    caloHitParameters.m_nCellInteractionLengths.Get()<<" thickness: "<<caloHitParameters.m_cellThickness.Get()<<std::endl;
}

//------------------------------------------------------------------------------------------------------------------------------------------

int DDCaloHitCreator::GetNLayersFromEdge(const EVENT::CalorimeterHit *const pCaloHit) const
{
    // Calo hit coordinate calculations
    const float barrelMaximumRadius(this->GetMaximumRadius(pCaloHit, m_settings.m_hCalBarrelOuterSymmetry, m_settings.m_hCalBarrelOuterPhi0));
    const float endCapMaximumRadius(this->GetMaximumRadius(pCaloHit, m_settings.m_hCalEndCapInnerSymmetryOrder, m_settings.m_hCalEndCapInnerPhiCoordinate));
    const float caloHitAbsZ(std::fabs(pCaloHit->getPosition()[2]));

    // Distance from radial outer
    float radialDistanceToEdge(std::numeric_limits<float>::max());

    if (caloHitAbsZ < m_settings.m_eCalBarrelOuterZ)
    {
        radialDistanceToEdge = (m_settings.m_hCalBarrelOuterR - barrelMaximumRadius) / m_hCalBarrelLayerThickness;
    }
    else
    {
        radialDistanceToEdge = (m_settings.m_hCalEndCapOuterR - endCapMaximumRadius) / m_hCalEndCapLayerThickness;
    }

    // Distance from rear of endcap outer
    float rearDistanceToEdge(std::numeric_limits<float>::max());

    if (caloHitAbsZ >= m_settings.m_eCalBarrelOuterZ)
    {
        rearDistanceToEdge = (m_settings.m_hCalEndCapOuterZ - caloHitAbsZ) / m_hCalEndCapLayerThickness;
    }
    else
    {
        const float rearDistance((m_settings.m_eCalBarrelOuterZ - caloHitAbsZ) / m_hCalBarrelLayerThickness);

        if (rearDistance < m_settings.m_layersFromEdgeMaxRearDistance)
        {
            const float overlapDistance((m_settings.m_hCalEndCapOuterR - endCapMaximumRadius) / m_hCalEndCapLayerThickness);
            rearDistanceToEdge = std::max(rearDistance, overlapDistance);
        }
    }

    return static_cast<int>(std::min(radialDistanceToEdge, rearDistanceToEdge));
}

//------------------------------------------------------------------------------------------------------------------------------------------

float DDCaloHitCreator::GetMaximumRadius(const EVENT::CalorimeterHit *const pCaloHit, const unsigned int symmetryOrder, const float phi0) const
{
    const float *pCaloHitPosition(pCaloHit->getPosition());

    if (symmetryOrder <= 2)
        return std::sqrt((pCaloHitPosition[0] * pCaloHitPosition[0]) + (pCaloHitPosition[1] * pCaloHitPosition[1]));

    float maximumRadius(0.f);
    const float twoPi(2.f * M_PI);

    for (unsigned int i = 0; i < symmetryOrder; ++i)
    {
        const float phi = phi0 + i * twoPi / static_cast<float>(symmetryOrder);
        float radius = pCaloHitPosition[0] * std::cos(phi) + pCaloHitPosition[1] * std::sin(phi);

        if (radius > maximumRadius)
            maximumRadius = radius;
    }

    return maximumRadius;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

DDCaloHitCreator::Settings::Settings()
  : m_eCalCaloHitCollections( StringVector() ),
    m_hCalCaloHitCollections( StringVector() ),
    m_lCalCaloHitCollections( StringVector() ),
    m_lHCalCaloHitCollections( StringVector() ),
    m_muonCaloHitCollections( StringVector() ),
    m_eCalToMip(1.f),
    m_hCalToMip(1.f),
    m_muonToMip(1.f),
    m_eCalMipThreshold(0.f),
    m_hCalMipThreshold(0.f),
    m_muonMipThreshold(0.f),
    m_eCalToEMGeV(1.f),
    m_eCalToHadGeVBarrel(1.f),
    m_eCalToHadGeVEndCap(1.f),
    m_hCalToEMGeV(1.f),
    m_hCalToHadGeV(1.f),
    m_muonDigitalHits(1),
    m_muonHitEnergy(0.5f),
    m_maxHCalHitHadronicEnergy(10000.f),
    m_nOuterSamplingLayers(3),
    m_layersFromEdgeMaxRearDistance(250.f),
    m_hCalEndCapInnerSymmetryOrder(4),
    m_hCalEndCapInnerPhiCoordinate(0.f),
    m_stripSplittingOn(0),
    m_useEcalScLayers(0),
    m_eCalSiToMip(1.f),
    m_eCalScToMip(1.f),
    m_eCalSiMipThreshold(0.f),
    m_eCalScMipThreshold(0.f),
    m_eCalSiToEMGeV(1.f),
    m_eCalScToEMGeV(1.f),
    m_eCalSiToHadGeVBarrel(1.f),
    m_eCalScToHadGeVBarrel(1.f),
    m_eCalSiToHadGeVEndCap(1.f),
    m_eCalScToHadGeVEndCap(1.f),
    m_eCalBarrelOuterZ(0.f),
    m_hCalBarrelOuterZ(0.f),
    m_muonBarrelOuterZ(0.f),
    m_coilOuterR(0.f),
    m_eCalBarrelInnerPhi0(0.f),
    m_eCalBarrelInnerSymmetry(0.f),
    m_hCalBarrelInnerPhi0(0.f),
    m_hCalBarrelInnerSymmetry(0.f),
    m_muonBarrelInnerPhi0(0.f),
    m_muonBarrelInnerSymmetry(0.f),
    m_hCalEndCapOuterR(0.f),
    m_hCalEndCapOuterZ(0.f),
    m_hCalBarrelOuterR(0.f),
    m_hCalBarrelOuterPhi0(0.f),
    m_hCalBarrelOuterSymmetry(0.f),
    m_eCalBarrelNormalVector({0.0, 0.0, 1.0}),
    m_hCalBarrelNormalVector({0.0, 0.0, 1.0}),
    m_muonBarrelNormalVector({0.0, 0.0, 1.0})
{
}
