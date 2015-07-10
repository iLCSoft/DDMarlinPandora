/**
 *  @file   MarlinPandora/src/DDGeometryCreator.cc
 * 
 *  @brief  Implementation of the geometry creator class.
 * 
 *  $Log: $
 */

#include "marlin/Global.h"
#include "marlin/Processor.h"

#include "DDGeometryCreator.h"
#include "DDPandoraPFANewProcessor.h"

#include "DD4hep/LCDD.h"
#include "DD4hep/DD4hepUnits.h"
#include "DDRec/DetectorData.h"

//Forward declarations. See DDPandoraPFANewProcessor.cc
DD4hep::DDRec::LayeredCalorimeterData * getExtension(std::string detectorName);

std::vector<double> getTrackingRegionExtent();
  

DDGeometryCreator::DDGeometryCreator(const Settings &settings, const pandora::Pandora *const pPandora) :
    m_settings(settings),
    m_pPandora(pPandora)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

DDGeometryCreator::~DDGeometryCreator()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode DDGeometryCreator::CreateGeometry() const
{
  
  DD4hep::Geometry::LCDD& lcdd = DD4hep::Geometry::LCDD::getInstance();
  
    try
    {
        SubDetectorTypeMap subDetectorTypeMap;
        this->SetMandatorySubDetectorParameters(subDetectorTypeMap);

        SubDetectorNameMap subDetectorNameMap;
//         this->SetAdditionalSubDetectorParameters(subDetectorNameMap);

        std::string detectorName = lcdd.world().name();
        
        if (std::string::npos != detectorName.find("ILD"))
            this->SetILDSpecificGeometry(subDetectorTypeMap, subDetectorNameMap);

        for (SubDetectorTypeMap::const_iterator iter = subDetectorTypeMap.begin(), iterEnd = subDetectorTypeMap.end(); iter != iterEnd; ++iter)
            PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::Geometry::SubDetector::Create(*m_pPandora, iter->second));

        for (SubDetectorNameMap::const_iterator iter = subDetectorNameMap.begin(), iterEnd = subDetectorNameMap.end(); iter != iterEnd; ++iter)
            PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::Geometry::SubDetector::Create(*m_pPandora, iter->second));
    }
    catch (std::exception &exception)
    {
        streamlog_out(ERROR) << "Failure in marlin pandora geometry creator, exception: " << exception.what() << std::endl;
        throw exception;
    }

    return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DDGeometryCreator::SetMandatorySubDetectorParameters(SubDetectorTypeMap &subDetectorTypeMap) const
{
    PandoraApi::Geometry::SubDetector::Parameters eCalBarrelParameters, eCalEndCapParameters, hCalBarrelParameters, hCalEndCapParameters,
        muonBarrelParameters, muonEndCapParameters;

    this->SetDefaultSubDetectorParameters(*const_cast<DD4hep::DDRec::LayeredCalorimeterData*>(getExtension(m_settings.m_ecalBarrelName)), "ECalBarrel", pandora::ECAL_BARREL, eCalBarrelParameters);
    this->SetDefaultSubDetectorParameters(*const_cast<DD4hep::DDRec::LayeredCalorimeterData*>(getExtension(m_settings.m_ecalEndcapName)), "ECalEndCap", pandora::ECAL_ENDCAP, eCalEndCapParameters);
    this->SetDefaultSubDetectorParameters(*const_cast<DD4hep::DDRec::LayeredCalorimeterData*>(getExtension(m_settings.m_hcalBarrelName)), "HCalBarrel", pandora::HCAL_BARREL, hCalBarrelParameters);
    this->SetDefaultSubDetectorParameters(*const_cast<DD4hep::DDRec::LayeredCalorimeterData*>(getExtension(m_settings.m_hcalEndcapName)), "HCalEndCap", pandora::HCAL_ENDCAP, hCalEndCapParameters);
    this->SetDefaultSubDetectorParameters(*const_cast<DD4hep::DDRec::LayeredCalorimeterData*>(getExtension(m_settings.m_muonBarrelName)), "MuonBarrel", pandora::MUON_BARREL, muonBarrelParameters);
    this->SetDefaultSubDetectorParameters(*const_cast<DD4hep::DDRec::LayeredCalorimeterData*>(getExtension(m_settings.m_muonEndcapName)), "MuonEndCap", pandora::MUON_ENDCAP, muonEndCapParameters);

    subDetectorTypeMap[pandora::ECAL_BARREL] = eCalBarrelParameters;
    subDetectorTypeMap[pandora::ECAL_ENDCAP] = eCalEndCapParameters;
    subDetectorTypeMap[pandora::HCAL_BARREL] = hCalBarrelParameters;
    subDetectorTypeMap[pandora::HCAL_ENDCAP] = hCalEndCapParameters;
    subDetectorTypeMap[pandora::MUON_BARREL] = muonBarrelParameters;
    subDetectorTypeMap[pandora::MUON_ENDCAP] = muonEndCapParameters;

    PandoraApi::Geometry::SubDetector::Parameters trackerParameters;
    std::cout<<"BIG WARNING IN DDGeometryCreator.cc: DON'T MASK! SHOULD CHANGE TO OBTAIN tpc/TRACKER PARAMS FROM LCDD"<<std::endl;
//     const gear::TPCParameters &tpcParameters(marlin::Global::GEAR->getTPCParameters());
    trackerParameters.m_subDetectorName = "Tracker";
    trackerParameters.m_subDetectorType = pandora::INNER_TRACKER;
    trackerParameters.m_innerRCoordinate = getTrackingRegionExtent()[0];
    trackerParameters.m_innerZCoordinate = 0.f;
    trackerParameters.m_innerPhiCoordinate = 0.f;
    trackerParameters.m_innerSymmetryOrder = 0;
    trackerParameters.m_outerRCoordinate = getTrackingRegionExtent()[1];
    trackerParameters.m_outerZCoordinate = getTrackingRegionExtent()[3];
    trackerParameters.m_outerPhiCoordinate = 0.f;
    trackerParameters.m_outerSymmetryOrder = 0;
    trackerParameters.m_isMirroredInZ = true;
    trackerParameters.m_nLayers = 0;
    subDetectorTypeMap[pandora::INNER_TRACKER] = trackerParameters;

    PandoraApi::Geometry::SubDetector::Parameters coilParameters;
    
    ///FIXME:Implement a parameter for the Coil/Solenoid name
    ///NOTE: Is this the way to go, or should we go with reco structure?
    DD4hep::Geometry::LCDD & lcdd = DD4hep::Geometry::LCDD::getInstance();
    DD4hep::Geometry::DetElement coilDE = lcdd.detector(m_settings.m_coilName);
    //access the detelement and create a shape from the envelope since only minimal info needed
    DD4hep::Geometry::Tube coilTube = DD4hep::Geometry::Tube( coilDE.volume().solid() )  ;
    
    
    coilParameters.m_subDetectorName = "Coil";
    coilParameters.m_subDetectorType = pandora::COIL;
    coilParameters.m_innerRCoordinate = coilTube->GetRmin()/ dd4hep::mm; 
    coilParameters.m_innerZCoordinate = 0.f;
    coilParameters.m_innerPhiCoordinate = 0.f;
    coilParameters.m_innerSymmetryOrder = 0;
    coilParameters.m_outerRCoordinate = coilTube->GetRmax()/ dd4hep::mm;
    coilParameters.m_outerZCoordinate = coilTube->GetDZ()/ dd4hep::mm;
    coilParameters.m_outerPhiCoordinate = 0.f;
    coilParameters.m_outerSymmetryOrder = 0;
    coilParameters.m_isMirroredInZ = true;
    coilParameters.m_nLayers = 0;
    subDetectorTypeMap[pandora::COIL] = coilParameters;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DDGeometryCreator::SetAdditionalSubDetectorParameters(SubDetectorNameMap &subDetectorNameMap) const
{
  
  for (std::vector<std::string>::const_iterator iter = m_settings.m_ecalOtherNames.begin(), iterEnd = m_settings.m_ecalOtherNames.end();iter != iterEnd; ++iter){
  
    try
    {
        PandoraApi::Geometry::SubDetector::Parameters parameters;
        this->SetDefaultSubDetectorParameters(*const_cast<DD4hep::DDRec::LayeredCalorimeterData*>(getExtension(*iter)), *iter, pandora::SUB_DETECTOR_OTHER, parameters);
        subDetectorNameMap[parameters.m_subDetectorName.Get()] = parameters;
    }
    catch (std::runtime_error &exception)
    {
        streamlog_out(WARNING) << "Marlin pandora geometry creator during Other ECal construction: " << exception.what() << std::endl;
    }
  }
  
  for (std::vector<std::string>::const_iterator iter = m_settings.m_hcalOtherNames.begin(), iterEnd = m_settings.m_hcalOtherNames.end();iter != iterEnd; ++iter){
    
    try
    {
      PandoraApi::Geometry::SubDetector::Parameters parameters;
      this->SetDefaultSubDetectorParameters(*const_cast<DD4hep::DDRec::LayeredCalorimeterData*>(getExtension(*iter)), *iter, pandora::SUB_DETECTOR_OTHER, parameters);
      subDetectorNameMap[parameters.m_subDetectorName.Get()] = parameters;
    }
    catch (std::runtime_error &exception)
    {
      streamlog_out(WARNING) << "Marlin pandora geometry creator during Other HCal construction: " << exception.what() << std::endl;
    }
  }
  
  for (std::vector<std::string>::const_iterator iter = m_settings.m_muonOtherNames.begin(), iterEnd = m_settings.m_muonOtherNames.end();iter != iterEnd; ++iter){
    
    try
    {
      PandoraApi::Geometry::SubDetector::Parameters parameters;
      this->SetDefaultSubDetectorParameters(*const_cast<DD4hep::DDRec::LayeredCalorimeterData*>(getExtension(*iter)), *iter, pandora::SUB_DETECTOR_OTHER, parameters);
      subDetectorNameMap[parameters.m_subDetectorName.Get()] = parameters;
    }
    catch (std::runtime_error &exception)
    {
      streamlog_out(WARNING) << "Marlin pandora geometry creator during Other Muon construction: " << exception.what() << std::endl;
    }
  }

}

//------------------------------------------------------------------------------------------------------------------------------------------

void DDGeometryCreator::SetDefaultSubDetectorParameters(const DD4hep::DDRec::LayeredCalorimeterData &inputParameters, const std::string &subDetectorName,
    const pandora::SubDetectorType subDetectorType, PandoraApi::Geometry::SubDetector::Parameters &parameters) const
{
  const std::vector<DD4hep::DDRec::LayeredCalorimeterStruct::Layer>& layers= inputParameters.layers;

    parameters.m_subDetectorName = subDetectorName;
    parameters.m_subDetectorType = subDetectorType;
    parameters.m_innerRCoordinate = inputParameters.extent[0]/dd4hep::mm;
    parameters.m_innerZCoordinate = inputParameters.extent[2]/dd4hep::mm;
    parameters.m_innerPhiCoordinate = inputParameters.phi0/dd4hep::rad;
    parameters.m_innerSymmetryOrder = inputParameters.inner_symmetry;
    parameters.m_outerRCoordinate = inputParameters.extent[1]/dd4hep::mm;
    parameters.m_outerZCoordinate = inputParameters.extent[3]/dd4hep::mm;
    parameters.m_outerPhiCoordinate = inputParameters.phi0/dd4hep::rad;
    parameters.m_outerSymmetryOrder = inputParameters.outer_symmetry;
    parameters.m_isMirroredInZ = true;
    parameters.m_nLayers = layers.size();

    // ATTN Not always going to be correct for any optional subdetectors, but impact of this is negligible for ILD
    ///FIXME: Should get from geometry?
    const float radiationLength(((pandora::ECAL_BARREL == subDetectorType) || (pandora::ECAL_ENDCAP == subDetectorType)) ? m_settings.m_absorberRadLengthECal :
        ((pandora::HCAL_BARREL == subDetectorType) || (pandora::HCAL_ENDCAP == subDetectorType)) ? m_settings.m_absorberRadLengthHCal : m_settings.m_absorberRadLengthOther);
    const float interactionLength(((pandora::ECAL_BARREL == subDetectorType) || (pandora::ECAL_ENDCAP == subDetectorType)) ? m_settings.m_absorberIntLengthECal :
        ((pandora::HCAL_BARREL == subDetectorType) || (pandora::HCAL_ENDCAP == subDetectorType)) ? m_settings.m_absorberIntLengthHCal : m_settings.m_absorberIntLengthOther);

    for (int i = 0; i < layers.size(); ++i)
    {
        const DD4hep::DDRec::LayeredCalorimeterStruct::Layer & theLayer = layers.at(i);
        
        PandoraApi::Geometry::LayerParameters layerParameters;
        layerParameters.m_closestDistanceToIp = theLayer.distance/dd4hep::mm; //FIXME! IS THIS NEEDED? + (0.5 * (theLayer.thickness/dd4hep::mm + theLayer.absorberThickness/dd4hep::mm));
        layerParameters.m_nRadiationLengths = radiationLength * theLayer.absorberThickness/dd4hep::mm;
        layerParameters.m_nInteractionLengths = interactionLength * theLayer.absorberThickness/dd4hep::mm;
        parameters.m_layerParametersList.push_back(layerParameters);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode DDGeometryCreator::SetILDSpecificGeometry(SubDetectorTypeMap &subDetectorTypeMap, SubDetectorNameMap &subDetectorNameMap) const
{
    
  
    // Set positions of gaps in ILD detector and add information missing from GEAR parameters file
    try
    {
        const DD4hep::DDRec::LayeredCalorimeterData * hCalBarrelParameters = getExtension(m_settings.m_hcalBarrelName);
        subDetectorTypeMap[pandora::HCAL_BARREL].m_outerPhiCoordinate = hCalBarrelParameters->phi0/dd4hep::rad;
//CHECK! WAS         hCalBarrelParameters.getIntVal("Hcal_outer_polygon_phi0"); 
        subDetectorTypeMap[pandora::HCAL_BARREL].m_outerSymmetryOrder =  hCalBarrelParameters->outer_symmetry;
//CHECK! was         hCalBarrelParameters.getIntVal("Hcal_outer_polygon_order");
    }
    catch (std::runtime_error &)
    {
        // aLaVideauGeometry
        return this->SetILD_SDHCALSpecificGeometry(subDetectorTypeMap);
    }

    subDetectorTypeMap[pandora::ECAL_ENDCAP].m_innerSymmetryOrder = m_settings.m_eCalEndCapInnerSymmetryOrder;
    subDetectorTypeMap[pandora::ECAL_ENDCAP].m_innerPhiCoordinate = m_settings.m_eCalEndCapInnerPhiCoordinate;
    subDetectorTypeMap[pandora::ECAL_ENDCAP].m_outerSymmetryOrder = m_settings.m_eCalEndCapOuterSymmetryOrder;
    subDetectorTypeMap[pandora::ECAL_ENDCAP].m_outerPhiCoordinate = m_settings.m_eCalEndCapOuterPhiCoordinate;

    subDetectorTypeMap[pandora::HCAL_ENDCAP].m_innerSymmetryOrder = m_settings.m_hCalEndCapInnerSymmetryOrder;
    subDetectorTypeMap[pandora::HCAL_ENDCAP].m_innerPhiCoordinate = m_settings.m_hCalEndCapInnerPhiCoordinate;
    subDetectorTypeMap[pandora::HCAL_ENDCAP].m_outerSymmetryOrder = m_settings.m_hCalEndCapOuterSymmetryOrder;
    subDetectorTypeMap[pandora::HCAL_ENDCAP].m_outerPhiCoordinate = m_settings.m_hCalEndCapOuterPhiCoordinate;

    ///FIXME! THIS SHOULD BE CHANGED
    subDetectorNameMap["HCalRing"].m_innerSymmetryOrder = m_settings.m_hCalRingInnerSymmetryOrder;
    subDetectorNameMap["HCalRing"].m_innerPhiCoordinate = m_settings.m_hCalRingInnerPhiCoordinate;
    subDetectorNameMap["HCalRing"].m_outerSymmetryOrder = m_settings.m_hCalRingOuterSymmetryOrder;
    subDetectorNameMap["HCalRing"].m_outerPhiCoordinate = m_settings.m_hCalRingOuterPhiCoordinate;

    // Gaps in detector active material
    PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, this->CreateHCalBarrelBoxGaps());
    PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, this->CreateHCalEndCapBoxGaps());
    PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, this->CreateHCalBarrelConcentricGaps());

    return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode DDGeometryCreator::SetILD_SDHCALSpecificGeometry(SubDetectorTypeMap &subDetectorTypeMap) const
{
    // Non-default values (and those missing from GEAR parameters file)...
    // The following 2 parameters have no sense for Videau Geometry, set them to 0
    subDetectorTypeMap[pandora::HCAL_BARREL].m_outerPhiCoordinate = 0;
    subDetectorTypeMap[pandora::HCAL_BARREL].m_outerSymmetryOrder = 0;

    // Endcap is identical to standard ILD geometry, only HCAL barrel is different
    subDetectorTypeMap[pandora::ECAL_ENDCAP].m_innerSymmetryOrder = m_settings.m_eCalEndCapInnerSymmetryOrder;
    subDetectorTypeMap[pandora::ECAL_ENDCAP].m_innerPhiCoordinate = m_settings.m_eCalEndCapInnerPhiCoordinate;
    subDetectorTypeMap[pandora::ECAL_ENDCAP].m_outerSymmetryOrder = m_settings.m_eCalEndCapOuterSymmetryOrder;
    subDetectorTypeMap[pandora::ECAL_ENDCAP].m_outerPhiCoordinate = m_settings.m_eCalEndCapOuterPhiCoordinate;

    subDetectorTypeMap[pandora::HCAL_ENDCAP].m_innerSymmetryOrder = m_settings.m_hCalEndCapInnerSymmetryOrder;
    subDetectorTypeMap[pandora::HCAL_ENDCAP].m_innerPhiCoordinate = m_settings.m_hCalEndCapInnerPhiCoordinate;
    subDetectorTypeMap[pandora::HCAL_ENDCAP].m_outerSymmetryOrder = m_settings.m_hCalEndCapOuterSymmetryOrder;
    subDetectorTypeMap[pandora::HCAL_ENDCAP].m_outerPhiCoordinate = m_settings.m_hCalEndCapOuterPhiCoordinate;

    // TODO implement gaps between modules

    return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode DDGeometryCreator::CreateHCalBarrelBoxGaps() const
{
    const std::string detectorName("CLIC"); //FIXME
    
    
    
    const DD4hep::DDRec::LayeredCalorimeterData * hCalBarrelParameters = getExtension(m_settings.m_hcalBarrelName);

    const unsigned int innerSymmetryOrder(hCalBarrelParameters->inner_symmetry);
    const unsigned int outerSymmetryOrder(hCalBarrelParameters->outer_symmetry);

    if ((0 == innerSymmetryOrder) || (2 != outerSymmetryOrder / innerSymmetryOrder))
    {
        streamlog_out(ERROR) << " Detector " << detectorName << " doesn't conform to expected ILD-specific geometry" << std::endl;
        return pandora::STATUS_CODE_INVALID_PARAMETER;
    }

    const float innerRadius(hCalBarrelParameters->extent[0]/dd4hep::mm);
    const float outerRadius(hCalBarrelParameters->extent[1]/dd4hep::mm);
    const float outerZ(hCalBarrelParameters->extent[3]/dd4hep::mm);
    const float phi0(hCalBarrelParameters->phi0/dd4hep::rad);

    const float staveGap(hCalBarrelParameters->gap0/dd4hep::mm);
    PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, this->CreateRegularBoxGaps(innerSymmetryOrder, phi0, innerRadius, outerRadius,
        -outerZ, outerZ, staveGap));

    const float outerPseudoPhi0(M_PI / static_cast<float>(innerSymmetryOrder));
    const float cosOuterPseudoPhi0(std::cos(outerPseudoPhi0));

    if ((0 == outerPseudoPhi0) || (0.f == cosOuterPseudoPhi0))
    {
        streamlog_out(ERROR) << " Detector " << detectorName << " doesn't conform to expected ILD-specific geometry" << std::endl;
        return pandora::STATUS_CODE_INVALID_PARAMETER;
    }

    const float middleStaveGap(hCalBarrelParameters->gap1/dd4hep::mm);

    PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, this->CreateRegularBoxGaps(innerSymmetryOrder, outerPseudoPhi0,
        innerRadius / cosOuterPseudoPhi0, outerRadius, -outerZ, outerZ, middleStaveGap));

    return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode DDGeometryCreator::CreateHCalEndCapBoxGaps() const
{
  
    
    const DD4hep::DDRec::LayeredCalorimeterData * hCalEndCapParameters = getExtension(m_settings.m_hcalEndcapName);

    const float staveGap(hCalEndCapParameters->gap0/dd4hep::mm);
    const float innerRadius(hCalEndCapParameters->extent[0]/dd4hep::mm);
    const float outerRadius(hCalEndCapParameters->extent[1]/dd4hep::mm);
    const float innerZ(hCalEndCapParameters->extent[2]/dd4hep::mm);
    const float outerZ(hCalEndCapParameters->extent[3]/dd4hep::mm);

    PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, this->CreateRegularBoxGaps(m_settings.m_hCalEndCapInnerSymmetryOrder,
        m_settings.m_hCalEndCapInnerPhiCoordinate, innerRadius, outerRadius, innerZ, outerZ, staveGap,
        pandora::CartesianVector(-innerRadius, 0, 0)));

    PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, this->CreateRegularBoxGaps(m_settings.m_hCalEndCapInnerSymmetryOrder,
        m_settings.m_hCalEndCapInnerPhiCoordinate, innerRadius, outerRadius, -outerZ, -innerZ, staveGap,
        pandora::CartesianVector(innerRadius, 0, 0)));

    return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode DDGeometryCreator::CreateHCalBarrelConcentricGaps() const
{
    
    
    const DD4hep::DDRec::LayeredCalorimeterData *hCalBarrelParameters = getExtension(m_settings.m_hcalBarrelName);
    const float gapWidth(hCalBarrelParameters->gap0/dd4hep::mm);

    PandoraApi::Geometry::ConcentricGap::Parameters gapParameters;

    gapParameters.m_minZCoordinate = -0.5f * gapWidth;
    gapParameters.m_maxZCoordinate =  0.5f * gapWidth;
    gapParameters.m_innerRCoordinate = hCalBarrelParameters->extent[0]/dd4hep::mm;
    gapParameters.m_innerPhiCoordinate = hCalBarrelParameters->phi0/dd4hep::rad;
    gapParameters.m_innerSymmetryOrder = hCalBarrelParameters->inner_symmetry;
    gapParameters.m_outerRCoordinate = hCalBarrelParameters->extent[1]/dd4hep::mm;
    gapParameters.m_outerPhiCoordinate = hCalBarrelParameters->phi0/dd4hep::rad;
    gapParameters.m_outerSymmetryOrder = hCalBarrelParameters->outer_symmetry;

    PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::Geometry::ConcentricGap::Create(*m_pPandora, gapParameters));

    return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode DDGeometryCreator::CreateRegularBoxGaps(unsigned int symmetryOrder, float phi0, float innerRadius, float outerRadius,
    float minZ, float maxZ, float gapWidth, pandora::CartesianVector vertexOffset) const
{
    const pandora::CartesianVector basicGapVertex(pandora::CartesianVector(-0.5f * gapWidth, innerRadius, minZ) + vertexOffset);
    const pandora::CartesianVector basicSide1(gapWidth, 0, 0);
    const pandora::CartesianVector basicSide2(0, outerRadius - innerRadius, 0);
    const pandora::CartesianVector basicSide3(0, 0, maxZ - minZ);

    for (unsigned int i = 0; i < symmetryOrder; ++i)
    {
        const float phi = phi0 + (2. * M_PI * static_cast<float>(i) / static_cast<float>(symmetryOrder));
        const float sinPhi(std::sin(phi));
        const float cosPhi(std::cos(phi));

        PandoraApi::Geometry::BoxGap::Parameters gapParameters;

        gapParameters.m_vertex = pandora::CartesianVector(cosPhi * basicGapVertex.GetX() + sinPhi * basicGapVertex.GetY(),
            -sinPhi * basicGapVertex.GetX() + cosPhi * basicGapVertex.GetY(), basicGapVertex.GetZ());
        gapParameters.m_side1 = pandora::CartesianVector(cosPhi * basicSide1.GetX() + sinPhi * basicSide1.GetY(),
            -sinPhi * basicSide1.GetX() + cosPhi * basicSide1.GetY(), basicSide1.GetZ());
        gapParameters.m_side2 = pandora::CartesianVector(cosPhi * basicSide2.GetX() + sinPhi * basicSide2.GetY(),
            -sinPhi * basicSide2.GetX() + cosPhi * basicSide2.GetY(), basicSide2.GetZ());
        gapParameters.m_side3 = pandora::CartesianVector(cosPhi * basicSide3.GetX() + sinPhi * basicSide3.GetY(),
            -sinPhi * basicSide3.GetX() + cosPhi * basicSide3.GetY(), basicSide3.GetZ());

        PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::Geometry::BoxGap::Create(*m_pPandora, gapParameters));
    }

    return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

DDGeometryCreator::Settings::Settings() :
    m_absorberRadLengthECal(1.f),
    m_absorberIntLengthECal(1.f),
    m_absorberRadLengthHCal(1.f),
    m_absorberIntLengthHCal(1.f),
    m_absorberRadLengthOther(1.f),
    m_absorberIntLengthOther(1.f),
    m_eCalEndCapInnerSymmetryOrder(4),
    m_eCalEndCapInnerPhiCoordinate(0.f),
    m_eCalEndCapOuterSymmetryOrder(8),
    m_eCalEndCapOuterPhiCoordinate(0.f),
    m_hCalEndCapInnerSymmetryOrder(4),
    m_hCalEndCapInnerPhiCoordinate(0.f),
    m_hCalEndCapOuterSymmetryOrder(16),
    m_hCalEndCapOuterPhiCoordinate(0.f),
    m_hCalRingInnerSymmetryOrder(8),
    m_hCalRingInnerPhiCoordinate(0.f),
    m_hCalRingOuterSymmetryOrder(16),
    m_hCalRingOuterPhiCoordinate(0.f)
{
}
