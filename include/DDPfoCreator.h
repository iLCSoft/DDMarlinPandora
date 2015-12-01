/**
 *  @file   DDMarlinPandora/include/DDPfoCreator.h
 * 
 *  @brief  Header file for the pfo creator class.
 * 
 *  $Log: $
 */

#ifndef DDPFOCREATOR_H
#define DDPFOCREATOR_H 1

#include "EVENT/LCEvent.h"

#include "Api/PandoraApi.h"

/**
 *  @brief  DDPfoCreator class
 */
class DDPfoCreator
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

        std::string     m_clusterCollectionName;                ///< The name of the cluster output collection
        std::string     m_pfoCollectionName;                    ///< The name of the pfo output collection
    };

    /**
     *  @brief  Constructor
     * 
     *  @param  settings the creator settings
     *  @param  pPandora address of the relevant pandora instance
     */
     DDPfoCreator(const Settings &settings, const pandora::Pandora *const pPandora);

    /**
     *  @brief  Destructor
     */
     ~DDPfoCreator();

    /**
     *  @brief  Create particle flow objects
     * 
     *  @param  pLCEvent the lcio event
     */    
    pandora::StatusCode CreateParticleFlowObjects(EVENT::LCEvent *pLCEvent);

private:
    const Settings          m_settings;                         ///< The pfo creator settings
    const pandora::Pandora *m_pPandora;                         ///< Address of the pandora object from which to extract the pfos
};

#endif // #ifndef DDPFOCREATOR_H
