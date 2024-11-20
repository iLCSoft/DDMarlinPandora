/**
 *  @file   DDMarlinPandora/src/DDGeometryCreator.cc
 * 
 *  @brief  Implementation of the geometry creator class.
 * 
 *  $Log: $
 */

#include "marlin/Global.h"
#include "marlin/Processor.h"

#include "DDGeometryCreatorALLEGRO.h"

#include "DD4hep/Detector.h"
#include "DD4hep/DD4hepUnits.h"
#include "DDRec/DetectorData.h"
#include "DD4hep/DetType.h"
#include "DD4hep/DetectorSelector.h"


#include <utility>

//Forward declarations. See DDPandoraPFANewProcessor.cc
// dd4hep::rec::LayeredCalorimeterData * getExtension(std::string detectorName);
dd4hep::rec::LayeredCalorimeterData * getExtension(unsigned int includeFlag, unsigned int excludeFlag=0);

std::vector<double> getTrackingRegionExtent();
  

DDGeometryCreatorALLEGRO::DDGeometryCreatorALLEGRO(const Settings &settings, const pandora::Pandora *const pPandora) : DDGeometryCreator(settings,pPandora)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

DDGeometryCreatorALLEGRO::~DDGeometryCreatorALLEGRO()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode DDGeometryCreatorALLEGRO::CreateGeometry() const
{
    try
    {
        SubDetectorTypeMap subDetectorTypeMap;
        this->SetMandatorySubDetectorParameters(subDetectorTypeMap);

        streamlog_out(DEBUG) << "Creating geometry for ALLEGRO detector"<< std::endl;

        for (SubDetectorTypeMap::const_iterator iter = subDetectorTypeMap.begin(), iterEnd = subDetectorTypeMap.end(); iter != iterEnd; ++iter)
            PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::Geometry::SubDetector::Create(m_pPandora, iter->second));
    }
    catch (std::exception &exception)
    {
        streamlog_out(ERROR) << "Failure in marlin pandora geometry creator, exception: " << exception.what() << std::endl;
        throw exception;
    }

    return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DDGeometryCreatorALLEGRO::SetMandatorySubDetectorParameters(SubDetectorTypeMap &subDetectorTypeMap) const
{
    PandoraApi::Geometry::SubDetector::Parameters eCalBarrelParameters, eCalEndCapParameters, hCalBarrelParameters, hCalEndCapParameters,
        muonBarrelParameters, muonEndCapParameters;

    this->SetDefaultSubDetectorParameters(*const_cast<dd4hep::rec::LayeredCalorimeterData*>(getExtension( ( dd4hep::DetType::CALORIMETER | dd4hep::DetType::ELECTROMAGNETIC | dd4hep::DetType::BARREL), ( dd4hep::DetType::AUXILIARY  |  dd4hep::DetType::FORWARD ) )), "ECalBarrel", pandora::ECAL_BARREL, eCalBarrelParameters);
    this->SetDefaultSubDetectorParameters(*const_cast<dd4hep::rec::LayeredCalorimeterData*>(getExtension( ( dd4hep::DetType::CALORIMETER | dd4hep::DetType::ELECTROMAGNETIC | dd4hep::DetType::ENDCAP), ( dd4hep::DetType::AUXILIARY  |  dd4hep::DetType::FORWARD ) )), "ECalEndCap", pandora::ECAL_ENDCAP, eCalEndCapParameters);
    this->SetDefaultSubDetectorParameters(*const_cast<dd4hep::rec::LayeredCalorimeterData*>(getExtension( ( dd4hep::DetType::CALORIMETER | dd4hep::DetType::HADRONIC | dd4hep::DetType::BARREL), ( dd4hep::DetType::AUXILIARY  |  dd4hep::DetType::FORWARD ) )), "HCalBarrel", pandora::HCAL_BARREL, hCalBarrelParameters);
    this->SetDefaultSubDetectorParameters(*const_cast<dd4hep::rec::LayeredCalorimeterData*>(getExtension( ( dd4hep::DetType::CALORIMETER | dd4hep::DetType::HADRONIC | dd4hep::DetType::ENDCAP), ( dd4hep::DetType::AUXILIARY  |  dd4hep::DetType::FORWARD ) )), "HCalEndCap", pandora::HCAL_ENDCAP, hCalEndCapParameters);
    this->SetDefaultSubDetectorParameters(*const_cast<dd4hep::rec::LayeredCalorimeterData*>(getExtension( ( dd4hep::DetType::CALORIMETER | dd4hep::DetType::MUON | dd4hep::DetType::BARREL), ( dd4hep::DetType::AUXILIARY  |  dd4hep::DetType::FORWARD ) )), "MuonBarrel", pandora::MUON_BARREL, muonBarrelParameters);
    this->SetDefaultSubDetectorParameters(*const_cast<dd4hep::rec::LayeredCalorimeterData*>(getExtension( ( dd4hep::DetType::CALORIMETER | dd4hep::DetType::MUON | dd4hep::DetType::ENDCAP), ( dd4hep::DetType::AUXILIARY  |  dd4hep::DetType::FORWARD ) )), "MuonEndCap", pandora::MUON_ENDCAP, muonEndCapParameters);

    subDetectorTypeMap[pandora::ECAL_BARREL] = eCalBarrelParameters;
    subDetectorTypeMap[pandora::ECAL_ENDCAP] = eCalEndCapParameters;
    subDetectorTypeMap[pandora::HCAL_BARREL] = hCalBarrelParameters;
    subDetectorTypeMap[pandora::HCAL_ENDCAP] = hCalEndCapParameters;
    subDetectorTypeMap[pandora::MUON_BARREL] = muonBarrelParameters;
    subDetectorTypeMap[pandora::MUON_ENDCAP] = muonEndCapParameters;

/*
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
*/
    /*
    ///FIXME:Implement a parameter for the Coil/Solenoid name
    ///NOTE: Is this the way to go, or should we go with reco structure?
    try{
        PandoraApi::Geometry::SubDetector::Parameters coilParameters;

        const dd4hep::rec::LayeredCalorimeterData * coilExtension= getExtension( ( dd4hep::DetType::COIL ) );



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
    */
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DDGeometryCreatorALLEGRO::SetDefaultSubDetectorParameters(const dd4hep::rec::LayeredCalorimeterData &inputParameters, const std::string &subDetectorName,
    const pandora::SubDetectorType subDetectorType, PandoraApi::Geometry::SubDetector::Parameters &parameters) const
{
  const std::vector<dd4hep::rec::LayeredCalorimeterStruct::Layer>& layers= inputParameters.layers;

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
        const dd4hep::rec::LayeredCalorimeterStruct::Layer & theLayer = layers.at(i);
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

