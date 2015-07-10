/**
 *  @file   MarlinPandora/include/DDGeometryCreator.h
 * 
 *  @brief  Header file for the geometry creator class.
 * 
 *  $Log: $
 */

#ifndef DDGEOMETRY_CREATOR_H
#define DDGEOMETRY_CREATOR_H 1

#include "Api/PandoraApi.h"

#include "DDRec/DetectorData.h"

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  DDGeometryCreator class
 */
class DDGeometryCreator
{
public:
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

        float           m_absorberRadLengthECal;                ///< The absorber radiation length in the ECal
        float           m_absorberIntLengthECal;                ///< The absorber interaction length in the ECal
        float           m_absorberRadLengthHCal;                ///< The absorber radiation length in the HCal
        float           m_absorberIntLengthHCal;                ///< The absorber interaction length in the HCal
        float           m_absorberRadLengthOther;               ///< The absorber radiation length in other detector regions
        float           m_absorberIntLengthOther;               ///< The absorber interaction length in other detector regions

        int             m_eCalEndCapInnerSymmetryOrder;         ///< ECal end cap inner symmetry order (missing from ILD gear files)
        float           m_eCalEndCapInnerPhiCoordinate;         ///< ECal end cap inner phi coordinate (missing from ILD gear files)
        int             m_eCalEndCapOuterSymmetryOrder;         ///< ECal end cap outer symmetry order (missing from ILD gear files)
        float           m_eCalEndCapOuterPhiCoordinate;         ///< ECal end cap outer phi coordinate (missing from ILD gear files)

        int             m_hCalEndCapInnerSymmetryOrder;         ///< HCal end cap inner symmetry order (missing from ILD gear files)
        float           m_hCalEndCapInnerPhiCoordinate;         ///< HCal end cap inner phi coordinate (missing from ILD gear files)
        int             m_hCalEndCapOuterSymmetryOrder;         ///< HCal end cap outer symmetry order (missing from ILD gear files)
        float           m_hCalEndCapOuterPhiCoordinate;         ///< HCal end cap outer phi coordinate (missing from ILD gear files)

        int             m_hCalRingInnerSymmetryOrder;           ///< HCal ring inner symmetry order (missing from ILD gear files)
        float           m_hCalRingInnerPhiCoordinate;           ///< HCal ring inner phi coordinate (missing from ILD gear files)
        int             m_hCalRingOuterSymmetryOrder;           ///< HCal ring outer symmetry order (missing from ILD gear files)
        float           m_hCalRingOuterPhiCoordinate;           ///< HCal ring outer phi coordinate (missing from ILD gear files)
        
        ///Added by Nikiforos
        std::string     m_vertexBarrelDetectorName;               ///< The vertex barrel detector name 
        std::vector<std::string> m_barrelTrackerNames;      ///< The Barrel Tracking detector names 
        std::vector<std::string> m_endcapTrackerNames;      ///< The Endcap Tracking detector names 
        std::string     m_ecalBarrelName;                   ///< The ECal barrel detector name 
        std::string     m_ecalEndcapName;                   ///< The ECal endcap detector name 
        std::vector<std::string>  m_ecalOtherNames;         ///< Additional ECal detector names 
        
        std::string     m_hcalBarrelName;                   ///< The HCal barrel detector name 
        std::string     m_hcalEndcapName;                   ///< The HCal endcap detector name 
        std::vector<std::string>  m_hcalOtherNames;         ///< Additional HCal detector names 
        
        std::string     m_muonBarrelName;                   ///< The Muon barrel detector name 
        std::string     m_muonEndcapName;                   ///< The Muon endcap detector name 
        std::vector<std::string>  m_muonOtherNames;         ///< Additional Muon detector names 
        
        std::string     m_coilName;                   ///< The detector name for the coil
        
    };

    /**
     *  @brief  Constructor
     * 
     *  @param  settings the creator settings
     *  @param  pPandora address of the relevant pandora instance
     */
     DDGeometryCreator(const Settings &settings, const pandora::Pandora *const pPandora);

    /**
     *  @brief  Destructor
     */
     ~DDGeometryCreator();

    /**
     *  @brief  Create geometry
     */
    pandora::StatusCode CreateGeometry() const;

private:
    typedef std::map<pandora::SubDetectorType, PandoraApi::Geometry::SubDetector::Parameters> SubDetectorTypeMap;
    typedef std::map<std::string, PandoraApi::Geometry::SubDetector::Parameters> SubDetectorNameMap;

    /**
     *  @brief  Set mandatory sub detector parameters
     * 
     *  @param  subDetectorTypeMap the sub detector type map
     */
    void SetMandatorySubDetectorParameters(SubDetectorTypeMap &subDetectorTypeMap) const;

    /**
     *  @brief  Set additional sub detector parameters
     * 
     *  @param  subDetectorNameMap the sub detector name map (for smaller sub detectors, identified uniquely only by name)
     */
    void SetAdditionalSubDetectorParameters(SubDetectorNameMap &subDetectorNameMap) const;

    /**
     *  @brief  Set sub detector parameters to their gear default values
     * 
     *  @param  inputParameters input parameters, from gear
     *  @param  subDetectorName the sub detector name
     *  @param  subDetectorType the sub detector type
     *  @param  parameters the sub detector parameters
     */
    void SetDefaultSubDetectorParameters(const DD4hep::DDRec::LayeredCalorimeterData &inputParameters, const std::string &subDetectorName,
        const pandora::SubDetectorType subDetectorType, PandoraApi::Geometry::SubDetector::Parameters &parameters) const;

    /**
     *  @brief  Set positions of gaps in ILD detector and add information missing from GEAR parameters file
     * 
     *  @param  subDetectorTypeMap the sub detector type map
     *  @param  subDetectorNameMap the sub detector name map (for smaller sub detectors, identified uniquely only by name)
     */
    pandora::StatusCode SetILDSpecificGeometry(SubDetectorTypeMap &subDetectorTypeMap, SubDetectorNameMap &subDetectorNameMap) const;

    /**
     *  @brief  Add information missing from GEAR parameters file for ILD SDHCAL detector
     * 
     *  @param  subDetectorTypeMap the sub detector type map
     */
    pandora::StatusCode SetILD_SDHCALSpecificGeometry(SubDetectorTypeMap &subDetectorTypeMap) const;

    /**
     *  @brief  Specify positions of hcal barrel box gaps - ILD specific
     */
    pandora::StatusCode CreateHCalBarrelBoxGaps() const;

    /**
     *  @brief  Specify positions of hcal end cap box gaps - ILD specific
     */
    pandora::StatusCode CreateHCalEndCapBoxGaps() const;

    /**
     *  @brief  Specify positions of hcal barrel concentric polygon gaps - ILD specific
     */
    pandora::StatusCode CreateHCalBarrelConcentricGaps() const;

    /**
     *  @brief  Create box gaps at regular positions on polygonal prism, oriented along main z axis - ILD specific
     * 
     *  @param  symmetryOrder the pandora geometry parameters
     *  @param  phi0 the phi coordinate
     *  @param  innerRadius the inner r coordinate
     *  @param  outerRadius the outer r coordinate
     *  @param  minZ the minimum z coordinate
     *  @param  maxZ the maximum z coordinate
     *  @param  gapWidth the gap width
     *  @param  vertexOffset position offset for vertex that doesn't point back to origin of xy plane
     */
    pandora::StatusCode CreateRegularBoxGaps(unsigned int symmetryOrder, float phi0, float innerRadius, float outerRadius, float minZ,
        float maxZ, float gapWidth, pandora::CartesianVector vertexOffset = pandora::CartesianVector(0, 0, 0)) const;

    const Settings          m_settings;                     ///< The geometry creator settings
    const pandora::Pandora *m_pPandora;                     ///< Address of the pandora object to create the geometry
};

#endif // #ifndef GEOMETRY_CREATOR_H
