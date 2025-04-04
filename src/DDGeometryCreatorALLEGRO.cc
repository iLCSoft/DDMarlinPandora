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

    // FIXME! AD: currently ignoring tracker parameters since we are using truth tracks
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

    ///FIXME:Implement a parameter for the Coil/Solenoid name
    ///NOTE: Is this the way to go, or should we go with reco structure?
    try{
        PandoraApi::Geometry::SubDetector::Parameters coilParameters;

        /*
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
        */

        // FIXME! AD: below an ugly way is used to define the geometry of the Solenoid since we do not have the LayeredCalorimeterData for Coil in ALLEGRO
        //Get ECal Barrel extension
        const dd4hep::rec::LayeredCalorimeterData * eCalBarrelExtension= getExtension( ( dd4hep::DetType::CALORIMETER | dd4hep::DetType::ELECTROMAGNETIC | dd4hep::DetType::BARREL),
                                                                                     ( dd4hep::DetType::AUXILIARY  |  dd4hep::DetType::FORWARD ) );
        //Get HCal Barrel extension
        const dd4hep::rec::LayeredCalorimeterData * hCalBarrelExtension= getExtension( ( dd4hep::DetType::CALORIMETER | dd4hep::DetType::HADRONIC | dd4hep::DetType::BARREL),
                                                                                     ( dd4hep::DetType::AUXILIARY |  dd4hep::DetType::FORWARD ) );
        coilParameters.m_subDetectorName = "Coil";
        coilParameters.m_subDetectorType = pandora::COIL;
        coilParameters.m_innerRCoordinate = eCalBarrelExtension->extent[1]/ dd4hep::mm; // ECAL barrel outer R is assumed to be the inner R for the Solenoid
        coilParameters.m_innerZCoordinate = 0.f;
        coilParameters.m_innerPhiCoordinate = 0.f;
        coilParameters.m_innerSymmetryOrder = 0;
        coilParameters.m_outerRCoordinate = hCalBarrelExtension->extent[0]/ dd4hep::mm; // HCAL barrel inner R is assumed to be the outer R for the Solenoid
        coilParameters.m_outerZCoordinate = eCalBarrelExtension->extent[3]/ dd4hep::mm; // Solenoid outer Z is assumed to be the same as ECAL barrel outer Z
        coilParameters.m_outerPhiCoordinate = 0.f;
        coilParameters.m_outerSymmetryOrder = 0;
        coilParameters.m_isMirroredInZ = true;
        coilParameters.m_nLayers = 0;
        subDetectorTypeMap[pandora::COIL] = coilParameters;

    } catch ( std::exception & e ) {
        streamlog_out(ERROR) << "Failed to access COIL parameters: "<<e.what()<<std::endl;
    }
}
