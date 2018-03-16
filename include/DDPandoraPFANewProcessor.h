/**
 *  @file   DDMarlinPandora/include/DDPandoraPFANewProcessor.h
 * 
 *  @brief  Header file for the pandora pfa new processor class.
 * 
 *  $Log: $
 */

#ifndef DDPANDORAPFANEWPROCESSOR_H
#define DDPANDORAPFANEWPROCESSOR_H 1

#include "DDCaloHitCreator.h"
#include "DDGeometryCreator.h"
#include "DDMCParticleCreator.h"
#include "DDPfoCreator.h"
#include "DDTrackCreatorBase.h"

#include "marlin/Processor.h"

#include "DD4hep/DetType.h"
#include "DD4hep/DetectorSelector.h"


namespace pandora {class Pandora;}

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  DDPandoraPFANewProcessor class
 */
class DDPandoraPFANewProcessor : public marlin::Processor
{
public:
    typedef std::vector<float> FloatVector;
    typedef std::vector<std::string> StringVector;

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

        std::string     m_pandoraSettingsXmlFile = "";      ///< The pandora settings xml file

        float           m_innerBField = 0.0;                ///< The bfield in the main tracker, ecal and hcal, units Tesla
        float           m_muonBarrelBField = 0.0;           ///< The bfield in the muon barrel, units Tesla
        float           m_muonEndCapBField = 0.0;           ///< The bfield in the muon endcap, units Tesla
        bool            m_useDD4hepField = false;           ///< Whether to use the DD4hep field map instead of the values above

        FloatVector     m_inputEnergyCorrectionPoints{};    ///< The input energy points for non-linearity energy correction
        FloatVector     m_outputEnergyCorrectionPoints{};   ///< The output energy points for non-linearity energy correction
        
        // Software compensation parameters
        FloatVector     m_softCompParameters{};
        FloatVector     m_softCompEnergyDensityBins{};
        float           m_energyDensityFinalBin = 0.0;
        float           m_maxClusterEnergyToApplySoftComp = 100.0;
        float           m_minCleanHitEnergy = 0.5;
        float           m_minCleanHitEnergyFraction = 0.01;
        float           m_minCleanCorrectedHitEnergy = 0.1;
        
        
        ///ADDED BY NIKIFOROS
        //Detector names not needed anymore, accessed by det type flags
        std::string     m_trackCreatorName = "";            ///< The name of the DDTrackCreator implementation to use

        
    };

    /**
     *  @brief  Default constructor
     */
    DDPandoraPFANewProcessor();
    DDPandoraPFANewProcessor(const DDPandoraPFANewProcessor&) = delete;
    DDPandoraPFANewProcessor& operator=(const DDPandoraPFANewProcessor&) = delete;

    /**
     *  @brief  Create new processor
     */
    virtual Processor *newProcessor();

    /**
     *  @brief  Initialize, called at startup
     */
    virtual void init();

    /**
     *  @brief  Process run header
     *
     *  @param  pLCRunHeader the lc run header
     */
    virtual void processRunHeader(lcio::LCRunHeader *pLCRunHeader);

    /**
     *  @brief  Process event, main entry point
     *
     *  @param  pLCEvent the lc event
     */
    virtual void processEvent(EVENT::LCEvent *pLCEvent);

    /**
     *  @brief  Checks for event
     *
     *  @param  pLCEvent the lc event
     */
    virtual void check(EVENT::LCEvent *pLCEvent);

    /**
     *  @brief  End, called at shutdown
     */
    virtual void end();

    /**
     *  @brief  Get address of the pandora instance
     * 
     *  @return address of the pandora instance
     */
    const pandora::Pandora *GetPandora() const;

    /**
     *  @brief  Get address of the current lcio event
     * 
     *  @param  pPandora address of the relevant pandora instance
     * 
     *  @return address of the current lcio event
     */
    static const EVENT::LCEvent *GetCurrentEvent(const pandora::Pandora *const pPandora);

private:
    /**
     *  @brief  Register user algorithm factories, energy correction functions and particle id functions,
     *          insert user code here
     */
    pandora::StatusCode RegisterUserComponents() const;

    /**
     *  @brief  Process steering file parameters, insert user code here
     */
    void ProcessSteeringFile();

    /**
     *  @brief  Copy some steering parameters between settings objects
     */
    void FinaliseSteeringParameters();

    /**
     *  @brief  Reset the pandora pfa new processor
     */
    void Reset();

    pandora::Pandora                     *m_pPandora = NULL;                  ///< Address of the pandora instance
    DDCaloHitCreator                     *m_pCaloHitCreator = NULL;           ///< The calo hit creator
    DDGeometryCreator                    *m_pGeometryCreator = NULL;          ///< The geometry creator
    DDTrackCreatorBase                   *m_pTrackCreator = NULL;             ///< The track creator
    DDMCParticleCreator                  *m_pDDMCParticleCreator = NULL;      ///< The mc particle creator
    DDPfoCreator                         *m_pDDPfoCreator = NULL;             ///< The pfo creator

    Settings                              m_settings{};                       ///< The settings for the pandora pfa new processor
    DDCaloHitCreator::Settings            m_caloHitCreatorSettings{};         ///< The calo hit creator settings
    DDGeometryCreator::Settings           m_geometryCreatorSettings{};        ///< The geometry creator settings
    DDMCParticleCreator::Settings         m_mcParticleCreatorSettings{};      ///< The mc particle creator settings
    DDTrackCreatorBase::Settings          m_trackCreatorSettings{};           ///< The track creator settings
    DDPfoCreator::Settings                m_pfoCreatorSettings{};             ///< The pfo creator settings

    typedef std::map<const pandora::Pandora *, EVENT::LCEvent *> PandoraToLCEventMap;
    static PandoraToLCEventMap          m_pandoraToLCEventMap;              ///< The pandora to lc event map
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline marlin::Processor *DDPandoraPFANewProcessor::newProcessor()
{
    return new DDPandoraPFANewProcessor;
}

#endif // #ifndef DDPANDORAPFANEWPROCESSOR_H
