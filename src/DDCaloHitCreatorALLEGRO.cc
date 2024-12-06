/**
 *  @file   DDMarlinPandora/src/DDCaloHitCreator.cc
 * 
 *  @brief  Implementation of the calo hit creator class.
 * 
 *  $Log: $
 */

#include "DDCaloHitCreatorALLEGRO.h"

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


DDCaloHitCreatorALLEGRO::DDCaloHitCreatorALLEGRO(const Settings &settings, const pandora::Pandora *const pPandora) :
    DDCaloHitCreator(settings, pPandora)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

DDCaloHitCreatorALLEGRO::~DDCaloHitCreatorALLEGRO()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode DDCaloHitCreatorALLEGRO::CreateCaloHits(const EVENT::LCEvent *const pLCEvent)
{
    PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, this->CreateECalCaloHits(pLCEvent));
    PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, this->CreateHCalCaloHits(pLCEvent));
    PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, this->CreateMuonCaloHits(pLCEvent));

    return pandora::STATUS_CODE_SUCCESS;
}


//------------------------------------------------------------------------------------------------------------------------------------------

void DDCaloHitCreatorALLEGRO::GetCommonCaloHitProperties(const EVENT::CalorimeterHit *const pCaloHit, PandoraApi::CaloHit::Parameters &caloHitParameters) const
{
    const float *pCaloHitPosition(pCaloHit->getPosition());
    const pandora::CartesianVector positionVector(pCaloHitPosition[0], pCaloHitPosition[1], pCaloHitPosition[2]);

    // FIXME! AD: for ECAL the cell gemoetry should be pandora::POINTING with cellSize0 = DeltaEta and cellSize1 = DeltaPhi
    caloHitParameters.m_cellGeometry = pandora::RECTANGULAR;
    caloHitParameters.m_positionVector = positionVector;
    caloHitParameters.m_expectedDirection = positionVector.GetUnitVector();
    caloHitParameters.m_pParentAddress = (void*)pCaloHit;
    caloHitParameters.m_inputEnergy = pCaloHit->getEnergy();
    caloHitParameters.m_time = pCaloHit->getTime();
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DDCaloHitCreatorALLEGRO::GetEndCapCaloHitProperties(const EVENT::CalorimeterHit *const pCaloHit, const std::vector<dd4hep::rec::LayeredCalorimeterStruct::Layer> &layers,
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
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DDCaloHitCreatorALLEGRO::GetBarrelCaloHitProperties( const EVENT::CalorimeterHit *const pCaloHit,
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
}

//------------------------------------------------------------------------------------------------------------------------------------------

int DDCaloHitCreatorALLEGRO::GetNLayersFromEdge(const EVENT::CalorimeterHit *const pCaloHit) const
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

float DDCaloHitCreatorALLEGRO::GetMaximumRadius(const EVENT::CalorimeterHit *const pCaloHit, const unsigned int symmetryOrder, const float phi0) const
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
