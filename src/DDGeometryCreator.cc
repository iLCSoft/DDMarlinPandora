/**
 *  @file   DDMarlinPandora/src/DDGeometryCreator.cc
 * 
 *  @brief  Implementation of the geometry creator class.
 * 
 *  $Log: $
 */

#include "marlin/Global.h"
#include "marlin/Processor.h"

#include "DDGeometryCreator.h"

#include "DD4hep/LCDD.h"
#include "DD4hep/DD4hepUnits.h"
#include "DDRec/DetectorData.h"
#include "DD4hep/DetType.h"
#include "DD4hep/DetectorSelector.h"


#include <utility>

//Forward declarations. See DDPandoraPFANewProcessor.cc
// DD4hep::DDRec::LayeredCalorimeterData * getExtension(std::string detectorName);
DD4hep::DDRec::LayeredCalorimeterData * getExtension(unsigned int includeFlag, unsigned int excludeFlag=0);

std::vector<double> getTrackingRegionExtent();
  

DDGeometryCreator::DDGeometryCreator(const Settings &settings, const pandora::Pandora *const pPandora) :
    m_settings(settings),
    m_pPandora(*pPandora)
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
        this->SetAdditionalSubDetectorParameters(subDetectorNameMap);

        std::string detectorName = lcdd.header().name();

        streamlog_out(DEBUG) << "Creating geometry for detector " << detectorName<< std::endl;
        
        //Before it was checking the detector name for the "ILD" substring 
        if (m_settings.m_createGaps)
            this->SetILDSpecificGeometry(subDetectorTypeMap, subDetectorNameMap);

        for (SubDetectorTypeMap::const_iterator iter = subDetectorTypeMap.begin(), iterEnd = subDetectorTypeMap.end(); iter != iterEnd; ++iter)
            PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::Geometry::SubDetector::Create(m_pPandora, iter->second));

        for (SubDetectorNameMap::const_iterator iter = subDetectorNameMap.begin(), iterEnd = subDetectorNameMap.end(); iter != iterEnd; ++iter){
            streamlog_out(DEBUG) << "Creating geometry for additional subdetector " << iter->first<< std::endl;

            PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::Geometry::SubDetector::Create(m_pPandora, iter->second));
        }
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

    this->SetDefaultSubDetectorParameters(*const_cast<DD4hep::DDRec::LayeredCalorimeterData*>(getExtension( ( DD4hep::DetType::CALORIMETER | DD4hep::DetType::ELECTROMAGNETIC | DD4hep::DetType::BARREL), ( DD4hep::DetType::AUXILIARY  |  DD4hep::DetType::FORWARD ) )), "ECalBarrel", pandora::ECAL_BARREL, eCalBarrelParameters);
    this->SetDefaultSubDetectorParameters(*const_cast<DD4hep::DDRec::LayeredCalorimeterData*>(getExtension( ( DD4hep::DetType::CALORIMETER | DD4hep::DetType::ELECTROMAGNETIC | DD4hep::DetType::ENDCAP), ( DD4hep::DetType::AUXILIARY  |  DD4hep::DetType::FORWARD ) )), "ECalEndCap", pandora::ECAL_ENDCAP, eCalEndCapParameters);
    this->SetDefaultSubDetectorParameters(*const_cast<DD4hep::DDRec::LayeredCalorimeterData*>(getExtension( ( DD4hep::DetType::CALORIMETER | DD4hep::DetType::HADRONIC | DD4hep::DetType::BARREL), ( DD4hep::DetType::AUXILIARY  |  DD4hep::DetType::FORWARD ) )), "HCalBarrel", pandora::HCAL_BARREL, hCalBarrelParameters);
    this->SetDefaultSubDetectorParameters(*const_cast<DD4hep::DDRec::LayeredCalorimeterData*>(getExtension( ( DD4hep::DetType::CALORIMETER | DD4hep::DetType::HADRONIC | DD4hep::DetType::ENDCAP), ( DD4hep::DetType::AUXILIARY  |  DD4hep::DetType::FORWARD ) )), "HCalEndCap", pandora::HCAL_ENDCAP, hCalEndCapParameters);
    this->SetDefaultSubDetectorParameters(*const_cast<DD4hep::DDRec::LayeredCalorimeterData*>(getExtension( ( DD4hep::DetType::CALORIMETER | DD4hep::DetType::MUON | DD4hep::DetType::BARREL), ( DD4hep::DetType::AUXILIARY  |  DD4hep::DetType::FORWARD ) )), "MuonBarrel", pandora::MUON_BARREL, muonBarrelParameters);
    this->SetDefaultSubDetectorParameters(*const_cast<DD4hep::DDRec::LayeredCalorimeterData*>(getExtension( ( DD4hep::DetType::CALORIMETER | DD4hep::DetType::MUON | DD4hep::DetType::ENDCAP), ( DD4hep::DetType::AUXILIARY  |  DD4hep::DetType::FORWARD ) )), "MuonEndCap", pandora::MUON_ENDCAP, muonEndCapParameters);

    subDetectorTypeMap[pandora::ECAL_BARREL] = eCalBarrelParameters;
    subDetectorTypeMap[pandora::ECAL_ENDCAP] = eCalEndCapParameters;
    subDetectorTypeMap[pandora::HCAL_BARREL] = hCalBarrelParameters;
    subDetectorTypeMap[pandora::HCAL_ENDCAP] = hCalEndCapParameters;
    subDetectorTypeMap[pandora::MUON_BARREL] = muonBarrelParameters;
    subDetectorTypeMap[pandora::MUON_ENDCAP] = muonEndCapParameters;

    PandoraApi::Geometry::SubDetector::Parameters trackerParameters;

    trackerParameters.m_subDetectorName = "Tracker";
    trackerParameters.m_subDetectorType = pandora::INNER_TRACKER;
    trackerParameters.m_innerRCoordinate = getTrackingRegionExtent()[0];
    trackerParameters.m_innerZCoordinate = 0.f;
    trackerParameters.m_innerPhiCoordinate = 0.f;
    trackerParameters.m_innerSymmetryOrder = 0;
    trackerParameters.m_outerRCoordinate = getTrackingRegionExtent()[1];
    trackerParameters.m_outerZCoordinate = getTrackingRegionExtent()[2];
    trackerParameters.m_outerPhiCoordinate = 0.f;
    trackerParameters.m_outerSymmetryOrder = 0;
    trackerParameters.m_isMirroredInZ = true;
    trackerParameters.m_nLayers = 0;
    subDetectorTypeMap[pandora::INNER_TRACKER] = trackerParameters;

    
    ///FIXME:Implement a parameter for the Coil/Solenoid name
    ///NOTE: Is this the way to go, or should we go with reco structure?
    try{
        
        PandoraApi::Geometry::SubDetector::Parameters coilParameters;

        const DD4hep::DDRec::LayeredCalorimeterData * coilExtension= getExtension( ( DD4hep::DetType::COIL ) );   



        coilParameters.m_subDetectorName = "Coil";
        coilParameters.m_subDetectorType = pandora::COIL;
        coilParameters.m_innerRCoordinate = coilExtension->extent[0]/ dd4hep::mm; 
        coilParameters.m_innerZCoordinate = 0.f;
        coilParameters.m_innerPhiCoordinate = 0.f;
        coilParameters.m_innerSymmetryOrder = 0;
        coilParameters.m_outerRCoordinate = coilExtension->extent[1]/ dd4hep::mm;
        coilParameters.m_outerZCoordinate = coilExtension->extent[3]/ dd4hep::mm;
        coilParameters.m_outerPhiCoordinate = 0.f;
        coilParameters.m_outerSymmetryOrder = 0;
        coilParameters.m_isMirroredInZ = true;
        coilParameters.m_nLayers = 0;
        subDetectorTypeMap[pandora::COIL] = coilParameters;
    } catch ( std::exception & e ) {
        
        streamlog_out(ERROR) << "Failed to access COIL parameters: "<<e.what()<<std::endl;
        
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DDGeometryCreator::SetAdditionalSubDetectorParameters(SubDetectorNameMap &subDetectorNameMap) const
{
    
  DD4hep::Geometry::LCDD & lcdd = DD4hep::Geometry::LCDD::getInstance();
  const std::vector< DD4hep::Geometry::DetElement>& theECalOtherDetectors = DD4hep::Geometry::DetectorSelector(lcdd).detectors(DD4hep::DetType::CALORIMETER | DD4hep::DetType::ELECTROMAGNETIC | DD4hep::DetType::AUXILIARY ) ;
  for (std::vector< DD4hep::Geometry::DetElement>::const_iterator iter = theECalOtherDetectors.begin(), iterEnd = theECalOtherDetectors.end();iter != iterEnd; ++iter){
   try
    {
        const DD4hep::Geometry::DetElement& theDetector = *iter;
        const DD4hep::DDRec::LayeredCalorimeterData * theExtension = theDetector.extension<DD4hep::DDRec::LayeredCalorimeterData>();

        PandoraApi::Geometry::SubDetector::Parameters parameters;
        this->SetDefaultSubDetectorParameters(*const_cast<DD4hep::DDRec::LayeredCalorimeterData*>(theExtension),theDetector.name(), pandora::SUB_DETECTOR_OTHER, parameters);
        subDetectorNameMap[parameters.m_subDetectorName.Get()] = parameters;
    }
    catch (std::runtime_error &exception)
    {
        streamlog_out(WARNING) << "Marlin pandora geometry creator during Other ECal construction: " << exception.what() << std::endl;
    }
  }
  
  const std::vector< DD4hep::Geometry::DetElement>& theHCalOtherDetectors = DD4hep::Geometry::DetectorSelector(lcdd).detectors(DD4hep::DetType::CALORIMETER | DD4hep::DetType::HADRONIC | DD4hep::DetType::AUXILIARY ) ;
  for (std::vector< DD4hep::Geometry::DetElement>::const_iterator iter = theHCalOtherDetectors.begin(), iterEnd = theHCalOtherDetectors.end();iter != iterEnd; ++iter){
   try
    {
        const DD4hep::Geometry::DetElement& theDetector = *iter;
        const DD4hep::DDRec::LayeredCalorimeterData * theExtension = theDetector.extension<DD4hep::DDRec::LayeredCalorimeterData>();

        PandoraApi::Geometry::SubDetector::Parameters parameters;
        this->SetDefaultSubDetectorParameters(*const_cast<DD4hep::DDRec::LayeredCalorimeterData*>(theExtension),theDetector.name(), pandora::SUB_DETECTOR_OTHER, parameters);
        subDetectorNameMap[parameters.m_subDetectorName.Get()] = parameters;
    }
    catch (std::runtime_error &exception)
    {
        streamlog_out(WARNING) << "Marlin pandora geometry creator during Other HCal construction: " << exception.what() << std::endl;
    }
  }
  
  const std::vector< DD4hep::Geometry::DetElement>& theMuonOtherDetectors = DD4hep::Geometry::DetectorSelector(lcdd).detectors(DD4hep::DetType::CALORIMETER | DD4hep::DetType::MUON| DD4hep::DetType::AUXILIARY ) ;
  for (std::vector< DD4hep::Geometry::DetElement>::const_iterator iter = theMuonOtherDetectors.begin(), iterEnd = theMuonOtherDetectors.end();iter != iterEnd; ++iter){
   try
    {
        const DD4hep::Geometry::DetElement& theDetector = *iter;
        const DD4hep::DDRec::LayeredCalorimeterData * theExtension = theDetector.extension<DD4hep::DDRec::LayeredCalorimeterData>();

        PandoraApi::Geometry::SubDetector::Parameters parameters;
        this->SetDefaultSubDetectorParameters(*const_cast<DD4hep::DDRec::LayeredCalorimeterData*>(theExtension),theDetector.name(), pandora::SUB_DETECTOR_OTHER, parameters);
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
    parameters.m_innerPhiCoordinate = inputParameters.inner_phi0/dd4hep::rad;
    parameters.m_innerSymmetryOrder = inputParameters.inner_symmetry;
    parameters.m_outerRCoordinate = inputParameters.extent[1]/dd4hep::mm;
    parameters.m_outerZCoordinate = inputParameters.extent[3]/dd4hep::mm;
    parameters.m_outerPhiCoordinate = inputParameters.outer_phi0/dd4hep::rad;
    parameters.m_outerSymmetryOrder = inputParameters.outer_symmetry;
    parameters.m_isMirroredInZ = true;
    parameters.m_nLayers = layers.size();

    for (size_t i = 0; i< layers.size(); i++)
    {
        const DD4hep::DDRec::LayeredCalorimeterStruct::Layer & theLayer = layers.at(i);
        
        PandoraApi::Geometry::LayerParameters layerParameters;
        
        double totalNumberOfRadLengths = theLayer.inner_nRadiationLengths;
        double totalNumberOfIntLengths = theLayer.inner_nInteractionLengths;
        
        if(i>0){
            //Add the numbers from previous layer's outer side
            totalNumberOfRadLengths += layers.at(i-1).outer_nRadiationLengths;
            totalNumberOfIntLengths += layers.at(i-1).outer_nInteractionLengths;
            
        }
        
        layerParameters.m_closestDistanceToIp = (theLayer.distance+theLayer.inner_thickness)/dd4hep::mm; //Distance to center of sensitive element
        layerParameters.m_nRadiationLengths = totalNumberOfRadLengths;
        layerParameters.m_nInteractionLengths = totalNumberOfIntLengths;

        parameters.m_layerParametersVector.push_back(layerParameters);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode DDGeometryCreator::SetILDSpecificGeometry(SubDetectorTypeMap &/*subDetectorTypeMap*/, SubDetectorNameMap &/*subDetectorNameMap*/) const
{
    streamlog_out(DEBUG0) << " Building gaps in detector active material"<<std::endl;
    // Gaps in detector active material
    PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, this->CreateHCalBarrelBoxGaps());
    PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, this->CreateHCalEndCapBoxGaps());
    PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, this->CreateHCalBarrelConcentricGaps());

    return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------




pandora::StatusCode DDGeometryCreator::CreateHCalBarrelBoxGaps() const
{
    DD4hep::Geometry::LCDD& lcdd = DD4hep::Geometry::LCDD::getInstance();

    std::string detectorName = lcdd.header().name();
    
    const DD4hep::DDRec::LayeredCalorimeterData * hCalBarrelParameters = getExtension(( DD4hep::DetType::CALORIMETER | DD4hep::DetType::HADRONIC | DD4hep::DetType::BARREL), ( DD4hep::DetType::AUXILIARY  |  DD4hep::DetType::FORWARD ));

    const unsigned int innerSymmetryOrder(hCalBarrelParameters->inner_symmetry);
    const unsigned int outerSymmetryOrder(hCalBarrelParameters->outer_symmetry);

    if ((0 == innerSymmetryOrder) || (2 != outerSymmetryOrder / innerSymmetryOrder))
    {
        streamlog_out(ERROR) << " Detector " << detectorName << " doesn't conform to expected ILD-specific geometry: innerSymmetryOrder : "<<innerSymmetryOrder <<" outerSymmetryOrder: "<<outerSymmetryOrder<< std::endl;
        return pandora::STATUS_CODE_INVALID_PARAMETER;
    }

    const float innerRadius(hCalBarrelParameters->extent[0]/dd4hep::mm);
    const float outerRadius(hCalBarrelParameters->extent[1]/dd4hep::mm);
    const float outerZ(hCalBarrelParameters->extent[3]/dd4hep::mm);
    const float inner_phi0(hCalBarrelParameters->inner_phi0/dd4hep::rad);

    const float staveGap(hCalBarrelParameters->gap0/dd4hep::mm);
    PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, this->CreateRegularBoxGaps(innerSymmetryOrder, inner_phi0, innerRadius, outerRadius,
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
  
    
    const DD4hep::DDRec::LayeredCalorimeterData * hCalEndCapParameters = getExtension(( DD4hep::DetType::CALORIMETER | DD4hep::DetType::HADRONIC | DD4hep::DetType::ENDCAP), ( DD4hep::DetType::AUXILIARY  |  DD4hep::DetType::FORWARD ));

    const float staveGap(hCalEndCapParameters->gap0/dd4hep::mm);
    const float innerRadius(hCalEndCapParameters->extent[0]/dd4hep::mm);
    const float outerRadius(hCalEndCapParameters->extent[1]/dd4hep::mm);
    const float innerZ(hCalEndCapParameters->extent[2]/dd4hep::mm);
    const float outerZ(hCalEndCapParameters->extent[3]/dd4hep::mm);
    const unsigned int innerSymmetryOrder(hCalEndCapParameters->inner_symmetry);
    const float inner_phi0(hCalEndCapParameters->inner_phi0/dd4hep::rad);

    PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, this->CreateRegularBoxGaps(innerSymmetryOrder,
        inner_phi0, innerRadius, outerRadius, innerZ, outerZ, staveGap,
        pandora::CartesianVector(-innerRadius, 0, 0)));

    PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, this->CreateRegularBoxGaps(innerSymmetryOrder,
        inner_phi0, innerRadius, outerRadius, -outerZ, -innerZ, staveGap,
        pandora::CartesianVector(innerRadius, 0, 0)));

    return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode DDGeometryCreator::CreateHCalBarrelConcentricGaps() const
{
    
    
    const DD4hep::DDRec::LayeredCalorimeterData *hCalBarrelParameters = getExtension(( DD4hep::DetType::CALORIMETER | DD4hep::DetType::HADRONIC | DD4hep::DetType::BARREL), ( DD4hep::DetType::AUXILIARY  |  DD4hep::DetType::FORWARD ));
    const float gapWidth(hCalBarrelParameters->gap0/dd4hep::mm);

    PandoraApi::Geometry::ConcentricGap::Parameters gapParameters;

    gapParameters.m_minZCoordinate = -0.5f * gapWidth;
    gapParameters.m_maxZCoordinate =  0.5f * gapWidth;
    gapParameters.m_innerRCoordinate = hCalBarrelParameters->extent[0]/dd4hep::mm;
    gapParameters.m_innerPhiCoordinate = hCalBarrelParameters->inner_phi0/dd4hep::rad;
    gapParameters.m_innerSymmetryOrder = hCalBarrelParameters->inner_symmetry;
    gapParameters.m_outerRCoordinate = hCalBarrelParameters->extent[1]/dd4hep::mm;
    gapParameters.m_outerPhiCoordinate = hCalBarrelParameters->outer_phi0/dd4hep::rad;
    gapParameters.m_outerSymmetryOrder = hCalBarrelParameters->outer_symmetry;

    PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::Geometry::ConcentricGap::Create(m_pPandora, gapParameters));

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

        PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::Geometry::BoxGap::Create(m_pPandora, gapParameters));
    }

    return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

DDGeometryCreator::Settings::Settings() :
    m_createGaps(false)
{
}
