/**
 *  @file   DDMarlinPandora/include/DDCaloHitCreator.h
 * 
 *  @brief  Header file for the calo hit creator class.
 * 
 *  $Log: $
 */

#ifndef DDCALO_HIT_CREATOR_H
#define DDCALO_HIT_CREATOR_H 1

#include "EVENT/CalorimeterHit.h"
#include "EVENT/LCEvent.h"

#include <IMPL/CalorimeterHitImpl.h>

#include "Api/PandoraApi.h"

#include <DDRec/DetectorData.h>
#include <DD4hep/Detector.h>
#include <DD4hep/DetElement.h>



typedef std::vector<EVENT::CalorimeterHit *> CalorimeterHitVector;

/**
 *  @brief  DDCaloHitCreator class
 */
class DDCaloHitCreator
{
public:
    typedef std::vector<std::string> StringVector;
    typedef std::vector<float> FloatVector;

    /**
     *  @brief  Settings class
     */
    class Settings
    {
    public:
        /**
         *  @brief  Default constructor
         */
        Settings();

        StringVector    m_eCalCaloHitCollections;               ///< The ecal calorimeter hit collections
        StringVector    m_hCalCaloHitCollections;               ///< The hcal calorimeter hit collections
        StringVector    m_lCalCaloHitCollections;               ///< The lcal calorimeter hit collections
        StringVector    m_lHCalCaloHitCollections;              ///< The lhcal calorimeter hit collections
        StringVector    m_muonCaloHitCollections;               ///< The muon calorimeter hit collections

        //NN: Material properties variables removed; obtained from geometry
        
        float           m_eCalToMip;                            ///< The calibration from deposited ECal energy to mip
        float           m_hCalToMip;                            ///< The calibration from deposited HCal energy to mip
        float           m_muonToMip;                            ///< The calibration from deposited Muon energy to mip
        float           m_eCalMipThreshold;                     ///< Threshold for creating calo hits in the ECal, units mip
        float           m_hCalMipThreshold;                     ///< Threshold for creating calo hits in the HCal, units mip
        float           m_muonMipThreshold;                     ///< Threshold for creating calo hits in the HCal, units mip

        float           m_eCalToEMGeV;                          ///< The calibration from deposited ECal energy to EM energy
        float           m_eCalToHadGeVBarrel;                   ///< The calibration from deposited ECal barrel energy to hadronic energy
        float           m_eCalToHadGeVEndCap;                   ///< The calibration from deposited ECal endcap energy to hadronic energy
        float           m_hCalToEMGeV;                          ///< The calibration from deposited HCal energy to EM energy
        float           m_hCalToHadGeV;                         ///< The calibration from deposited HCal energy to hadronic energy
        int             m_muonDigitalHits;                      ///< Muon hits are treated as digital (energy from hit count)
        float           m_muonHitEnergy;                        ///< The energy for a digital muon calorimeter hit, units GeV

        float           m_maxHCalHitHadronicEnergy;             ///< The maximum hadronic energy allowed for a single hcal hit
        int             m_nOuterSamplingLayers;                 ///< Number of layers from edge for hit to be flagged as an outer layer hit
        float           m_layersFromEdgeMaxRearDistance;        ///< Maximum number of layers from candidate outer layer hit to rear of detector

        int             m_hCalEndCapInnerSymmetryOrder;         ///< HCal end cap inner symmetry order (missing from ILD00 gear file)
        float           m_hCalEndCapInnerPhiCoordinate;         ///< HCal end cap inner phi coordinate (missing from ILD00 gear file)

        // For Strip Splitting method and hybrid ECAL.
        int             m_stripSplittingOn;                     ///< To use SSA, this should be true (default is false)
        int             m_useEcalScLayers;                      ///< To use scintillator layers ~ hybrid ECAL, this should be true (default is false)
        float           m_eCalSiToMip;                          ///< The calibration from deposited Si-layer energy to mip
        float           m_eCalScToMip;                          ///< The calibration from deposited Sc-layer energy to mip
        float           m_eCalSiMipThreshold;                   ///< Threshold for creating calo hits in the Si-layers of ECAL, units mip
        float           m_eCalScMipThreshold;                   ///< Threshold for creating calo hits in the Sc-layers of ECAL, units mip
        float           m_eCalSiToEMGeV;                        ///< The calibration from deposited Si-layer energy to EM energy
        float           m_eCalScToEMGeV;                        ///< The calibration from deposited Sc-layer energy to EM energy
        float           m_eCalSiToHadGeVBarrel;                 ///< The calibration from deposited Si-layer energy on the enecaps to hadronic energy
        float           m_eCalScToHadGeVBarrel;                 ///< The calibration from deposited Sc-layer energy on the endcaps to hadronic energy
        float           m_eCalSiToHadGeVEndCap;                 ///< The calibration from deposited Si-layer energy on the enecaps to hadronic energy
        float           m_eCalScToHadGeVEndCap;                 ///< The calibration from deposited Sc-layer energy on the endcaps to hadronic energy
        
        
        ///ADDED BY NIKIFOROS
        //Note that names are not needed anymore since the detector elements can be accessed by type flags
        //Nikiforos: Moved from main class 
        
        float                         m_eCalBarrelOuterZ;                 ///< ECal barrel outer z coordinate
        float                         m_hCalBarrelOuterZ;                 ///< HCal barrel outer z coordinate
        float                         m_muonBarrelOuterZ;                 ///< Muon barrel outer z coordinate
        float                         m_coilOuterR;                       ///< Coil outer r coordinate
          
        float                         m_eCalBarrelInnerPhi0;              ///< ECal barrel inner phi0 coordinate
        unsigned int                  m_eCalBarrelInnerSymmetry;          ///< ECal barrel inner symmetry order
        float                         m_hCalBarrelInnerPhi0;              ///< HCal barrel inner phi0 coordinate
        unsigned int                  m_hCalBarrelInnerSymmetry;          ///< HCal barrel inner symmetry order
        float                         m_muonBarrelInnerPhi0;              ///< Muon barrel inner phi0 coordinate
        unsigned int                  m_muonBarrelInnerSymmetry;          ///< Muon barrel inner symmetry order
          
        float                         m_hCalEndCapOuterR;                 ///< HCal endcap outer r coordinate
        float                         m_hCalEndCapOuterZ;                 ///< HCal endcap outer z coordinate
        float                         m_hCalBarrelOuterR;                 ///< HCal barrel outer r coordinate
        float                         m_hCalBarrelOuterPhi0;              ///< HCal barrel outer phi0 coordinate
        unsigned int                  m_hCalBarrelOuterSymmetry;          ///< HCal barrel outer symmetry order

    public:
      FloatVector m_eCalBarrelNormalVector;
      FloatVector m_hCalBarrelNormalVector;
      FloatVector m_muonBarrelNormalVector;


    };

    /**
     *  @brief  Constructor
     * 
     *  @param  settings the creator settings
     *  @param  pPandora address of the relevant pandora instance
     */
     DDCaloHitCreator(const Settings &settings, const pandora::Pandora *const pPandora);

    /**
     *  @brief  Destructor
     */
     ~DDCaloHitCreator();

    /**
     *  @brief  Create calo hits
     * 
     *  @param  pLCEvent the lcio event
     */    
    pandora::StatusCode CreateCaloHits(const EVENT::LCEvent *const pLCEvent);

    /**
     *  @brief  Get the calorimeter hit vector
     * 
     *  @return The calorimeter hit vector
     */
    const CalorimeterHitVector &GetCalorimeterHitVector() const;

    /**
     *  @brief  Reset the calo hit creator
     */
    void Reset();

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
     *  @brief  Create lcal calo hits
     * 
     *  @param  pLCEvent the lcio event
     */    
    pandora::StatusCode CreateLCalCaloHits(const EVENT::LCEvent *const pLCEvent);

    /**
     *  @brief  Create lhcal calo hits
     * 
     *  @param  pLCEvent the lcio event
     */
    pandora::StatusCode CreateLHCalCaloHits(const EVENT::LCEvent *const pLCEvent);

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

    const Settings                      m_settings;                         ///< The calo hit creator settings

    const pandora::Pandora &            m_pandora;                          ///< Reference to the pandora object to create calo hits

    float                               m_hCalBarrelLayerThickness;         ///< HCal barrel layer thickness
    float                               m_hCalEndCapLayerThickness;         ///< HCal endcap layer thickness

    CalorimeterHitVector                m_calorimeterHitVector;             ///< The calorimeter hit vector

    dd4hep::VolumeManager m_volumeManager; ///< DD4hep volume manager

};

//------------------------------------------------------------------------------------------------------------------------------------------

inline const CalorimeterHitVector &DDCaloHitCreator::GetCalorimeterHitVector() const
{
    return m_calorimeterHitVector;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void DDCaloHitCreator::Reset()
{
    m_calorimeterHitVector.clear();
}

#endif // #ifndef CALO_HIT_CREATOR_H
