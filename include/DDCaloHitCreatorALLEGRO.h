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
     *  @brief  Create ecal calo hits
     * 
     *  @param  pLCEvent the lcio event
     */
    pandora::StatusCode CreateECalCaloHits(const EVENT::LCEvent *const pLCEvent);

    /**
     *  @brief  Create hcal calo hits
     * 
     *  @param  pLCEvent the lcio event
     */
    pandora::StatusCode CreateHCalCaloHits(const EVENT::LCEvent *const pLCEvent);

    /**
     *  @brief  Create muon calo hits
     * 
     *  @param  pLCEvent the lcio event
     */
    pandora::StatusCode CreateMuonCaloHits(const EVENT::LCEvent *const pLCEvent);

    /**
     *  @brief  Get common calo hit properties: position, parent address, input energy and time
     * 
     *  @param  pCaloHit the lcio calorimeter hit
     *  @param  caloHitParameters the calo hit parameters to populate
     */
    void GetCommonCaloHitProperties(const EVENT::CalorimeterHit *const pCaloHit, PandoraApi::CaloHit::Parameters &caloHitParameters) const;

    /**
     *  @brief  Get end cap specific calo hit properties: cell size, absorber radiation and interaction lengths, normal vector
     * 
     *  @param  pCaloHit the lcio calorimeter hit
     *  @param  layers the vector of layers from DDRec extensions
     *  @param  caloHitParameters the calo hit parameters to populate
     *  @param  absorberCorrection to receive the absorber thickness correction for the mip equivalent energy
     */
    void GetEndCapCaloHitProperties(const EVENT::CalorimeterHit *const pCaloHit,  const std::vector<dd4hep::rec::LayeredCalorimeterStruct::Layer> &layers,
        PandoraApi::CaloHit::Parameters &caloHitParameters, float &absorberCorrection) const;

    /**
     *  @brief  Get barrel specific calo hit properties: cell size, absorber radiation and interaction lengths, normal vector
     * 
     *  @param  pCaloHit the lcio calorimeter hit
     *  @param  layers the vector of layers from DDRec extensions
     *  @param  barrelSymmetryOrder the barrel order of symmetry
     *  @param  caloHitParameters the calo hit parameters to populate
     *  @param  normalVector is the normalVector to the sensitive layers in local coordinates
     *  @param  absorberCorrection to receive the absorber thickness correction for the mip equivalent energy
     */
    void GetBarrelCaloHitProperties( const EVENT::CalorimeterHit *const pCaloHit,
				     const std::vector<dd4hep::rec::LayeredCalorimeterStruct::Layer> &layers,
				     unsigned int barrelSymmetryOrder,
				     PandoraApi::CaloHit::Parameters &caloHitParameters,
				     FloatVector const& normalVector,
				     float &absorberCorrection ) const;

    /**
     *  @brief  Get number of active layers from position of a calo hit to the edge of the detector
     * 
     *  @param  pCaloHit the lcio calorimeter hit
     */
    int GetNLayersFromEdge(const EVENT::CalorimeterHit *const pCaloHit) const;

    /**
     *  @brief  Get the maximum radius of a calo hit in a polygonal detector structure
     * 
     *  @param  pCaloHit the lcio calorimeter hit
     *  @param  symmetryOrder the symmetry order
     *  @param  phi0 the angular orientation
     * 
     *  @return the maximum radius
     */
    float GetMaximumRadius(const EVENT::CalorimeterHit *const pCaloHit, const unsigned int symmetryOrder, const float phi0) const;

};
#endif // #ifndef CALO_HIT_CREATOR_H
