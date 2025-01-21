/**
 *  @file   DDMarlinPandora/include/DDCaloHitCreatorALLEGRO.h
 * 
 *  @brief  Header file for the calo hit creator class.
 * 
 *  $Log: $
 */

#ifndef DDCALO_HIT_CREATORALLEGRO_H
#define DDCALO_HIT_CREATORALLEGRO_H 1

#include "EVENT/CalorimeterHit.h"
#include "EVENT/LCEvent.h"

#include <IMPL/CalorimeterHitImpl.h>

#include "Api/PandoraApi.h"

#include <DDRec/DetectorData.h>
#include <DD4hep/Detector.h>
#include <DD4hep/DetElement.h>

#include "DDCaloHitCreator.h"


/**
 *  @brief  DDCaloHitCreator class
 */
class DDCaloHitCreatorALLEGRO : public DDCaloHitCreator
{
public:
    typedef std::vector<std::string> StringVector;
    typedef std::vector<float> FloatVector;

    /**
     *  @brief  Constructor
     * 
     *  @param  settings the creator settings
     *  @param  pPandora address of the relevant pandora instance
     */
     DDCaloHitCreatorALLEGRO(const Settings &settings, const pandora::Pandora *const pPandora);

    /**
     *  @brief  Destructor
     */
     ~DDCaloHitCreatorALLEGRO();

    /**
     *  @brief  Create calo hits
     * 
     *  @param  pLCEvent the lcio event
     */    
    pandora::StatusCode CreateCaloHits(const EVENT::LCEvent *const pLCEvent);

private:
    /**
     *  @brief  Get common calo hit properties: position, parent address, input energy and time
     * 
     *  @param  pCaloHit the lcio calorimeter hit
     *  @param  caloHitParameters the calo hit parameters to populate
     */
    void GetCommonCaloHitProperties(const EVENT::CalorimeterHit *const pCaloHit, PandoraApi::CaloHit::Parameters &caloHitParameters) const;

};
#endif // #ifndef CALO_HIT_CREATOR_H
