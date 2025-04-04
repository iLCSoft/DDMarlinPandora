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
