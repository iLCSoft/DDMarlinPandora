/**
 *  @file   MarlinPandora/include/MCParticleCreator.h
 * 
 *  @brief  Header file for the mc particle creator class.
 * 
 *  $Log: $
 */

#ifndef MC_PARTICLE_CREATOR_H
#define MC_PARTICLE_CREATOR_H 1

#include "EVENT/LCEvent.h"

#include "Api/PandoraApi.h"

#include "DDCaloHitCreator.h"
#include "DDTrackCreator.h"

/**
 *  @brief  MCParticleCreator class
 */
class MCParticleCreator
{
public:
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

        StringVector    m_mcParticleCollections;                ///< The mc particle collections
        StringVector    m_lcCaloHitRelationCollections;         ///< The SimCaloHit to CaloHit particle relations
        StringVector    m_lcTrackRelationCollections;           ///< The SimTrackerHit to TrackerHit particle relations
    };

    /**
     *  @brief  Constructor
     * 
     *  @param  settings the creator settings
     *  @param  pPandora address of the relevant pandora instance
     */
     MCParticleCreator(const Settings &settings, const pandora::Pandora *const pPandora);

    /**
     *  @brief  Destructor
     */
     ~MCParticleCreator();

    /**
     *  @brief  Create MCParticles
     * 
     *  @param  pLCEvent the lcio event
     */    
    pandora::StatusCode CreateMCParticles(const EVENT::LCEvent *const pLCEvent) const;

    /**
     *  @brief  Create Track to mc particle relationships
     *
     *  @param  pLCEvent the lcio event
     *  @param  trackVector the vector containing all tracks successfully passed to pandora
     */
    pandora::StatusCode CreateTrackToMCParticleRelationships(const EVENT::LCEvent *const pLCEvent, const TrackVector &trackVector) const;

    /**
     *  @brief  Create calo hit to mc particle relationships
     *
     *  @param  pLCEvent the lcio event
     *  @param  calorimeterHitVector the vector containing all calorimeter hits successfully passed to pandora
     */
    pandora::StatusCode CreateCaloHitToMCParticleRelationships(const EVENT::LCEvent *const pLCEvent, const CalorimeterHitVector &calorimeterHitVector) const;

private:
    const Settings          m_settings;                         ///< The mc particle creator settings
    const pandora::Pandora *m_pPandora;                         ///< Address of the pandora object to create the mc particles
    const float             m_bField;                           ///< The bfield
};

#endif // #ifndef MC_PARTICLE_CREATOR_H
