/**
 *  @file   DDMarlinPandora/include/DDTrackCreatorALLEGRO.h
 *
 *  @brief  Header file for the ILD implementation of the track creator class.
 *
 *  $Log: $
 */

#ifndef DDTRACK_CREATOR_ALLEGRO_H
#define DDTRACK_CREATOR_ALLEGRO_H 1

#include "DDTrackCreatorBase.h"

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  DDTrackCreatorALLEGRO class
 */
class DDTrackCreatorALLEGRO : public DDTrackCreatorBase
{
public:

    /**
     *  @brief  Constructor
     *
     *  @param  settings the creator settings
     *  @param  pPandora address of the relevant pandora instance
     */
     DDTrackCreatorALLEGRO(const Settings &settings, const pandora::Pandora *const pPandora);

    /**
     *  @brief  Destructor
     */
     virtual ~DDTrackCreatorALLEGRO();

    /**
     *  @brief  Create tracks implementation, insert user code here. Detector specific
     *
     *  @param  pLCEvent the lcio event
     */
     pandora::StatusCode CreateTracks(EVENT::LCEvent *pLCEvent);



protected:


    // Detector-specific configuration variables

    // drift chamber
    float                   m_cosDch;                       ///< Cos(theta) value at end of dch
    float                   m_dchInnerR;                    ///< The dch inner radius
    float                   m_dchOuterR;                    ///< The dch outer radius
    float                   m_dchOuterZ;                    ///< The dch maximum z coordinate
    unsigned int            m_dchNLayers;                   ///< The number of DCH layers

    // wrapper

    float                   m_wrapperBarrelInnerR;          ///< Inner radius of the barrel wrapper
    float                   m_wrapperBarrelOuterR;          ///< Outer radius of the barrel wrapper
    float                   m_wrapperBarrelOuterZ;          ///< Maximum z of the barrel wrapper
    float                   m_wrapperEndcapInnerR;          ///< Inner radius of the endcap wrapper
    float                   m_wrapperEndcapOuterR;          ///< Outer radius of the endcap wrapper
    float                   m_wrapperEndcapInnerZ;          ///< Minimum |z| of the endcap wrapper
    float                   m_wrapperEndcapOuterZ;          ///< Maximum |z| of the endcap wrapper
    unsigned int            m_wrapperBarrelNLayers;         ///< Number of layers of Si wrapper in barrel
    unsigned int            m_wrapperEndcapNLayers;         ///< Number of layers of Si wrapper in endcap
    float m_cosWrapper;                                     ///< Cos(theta) for tracks pointing from IP to tracker edge
    // DoubleVector            m_barrelWrapperRPositions;      ///< List of barrel wrapper r positions
    // DoubleVector            m_barrelWrapperOuterZ;          ///< List of barrel wrapper outer |z|
    // DoubleVector            m_endcapWrapperInnerR;          ///< List of endcap wrapper inner radii
    // DoubleVector            m_endcapWrapperOuterR;          ///< List of endcap wrapper outer radii
    // DoubleVector            m_endcapWrapperZPositions;      ///< List of endcap wrapper z positions

    /**
     *  @brief  Whether track passes the quality cuts required in order to be used to form a pfo. Detector specific
     *
     *  @param  pTrack the lcio track
     *  @param  trackParameters the track parameters
     *
     *  @return boolean
     */

    virtual bool PassesQualityCuts(const EVENT::Track *const pTrack, const PandoraApi::Track::Parameters &trackParameters) const;

    /**
     *  @brief  Get number of hits in vertex of a track
     *
     *  @param  pTrack the lcio track
     *
     *  @return number of hits in vertex of a track
     */
    int GetNVertexHits(const EVENT::Track *const pTrack) const;

    /**
     *  @brief  Get number of hits in DCH of a track
     *
     *  @param  pTrack the lcio track
     *
     *  @return number of hits in DCH of a track
     */
    int GetNDchHits(const EVENT::Track *const pTrack) const;

    /**
     *  @brief  Get number of hits in Si wrapper of a track
     *
     *  @param  pTrack the lcio track
     *
     *  @return number of hits in Si wrapper of a track
     */
    int GetNSiWrapperHits(const EVENT::Track *const pTrack) const;

    /**
     *  @brief  Decide whether track reaches the ecal surface. Detector specific
     *
     *  @param  pTrack the lcio track
     *  @param  trackParameters the track parameters
     */
    virtual void TrackReachesECAL(const EVENT::Track *const pTrack, PandoraApi::Track::Parameters &trackParameters) const;

    /**
     *  @brief  Determine whether a track can be used to form a pfo under the following conditions:
     *          1) if the track proves to be associated with a cluster, OR
     *          2) if the track proves to have no cluster associations
     *          Detector specific
     *
     *  @param  pTrack the lcio track
     *  @param  trackParameters the track parameters
     */
    virtual void DefineTrackPfoUsage(const EVENT::Track *const pTrack, PandoraApi::Track::Parameters &trackParameters) const;

    /**
     *  @brief  Copy track states stored in lcio tracks to pandora track parameters
     *
     *  @param  pTrack the lcio track
     *  @param  trackParameters the track parameters
     */
    // inherited from base class
    // void GetTrackStates(const EVENT::Track *const pTrack, PandoraApi::Track::Parameters &trackParameters) const;

    /**
     *  @brief  Obtain track time when it reaches ECAL
     *
     *  @param  pTrack the lcio track
     */
    // inherited from base class
    // float CalculateTrackTimeAtCalorimeter(const EVENT::Track *const pTrack) const;

    /**
     *  @brief  Calculate possible second track state at the ECal Endcap
     *
     *  @param track lcio track
     *  @param trackParameters pandora LCTrackParameters
     */
    virtual void GetTrackStatesAtCalo( EVENT::Track *track, lc_content::LCTrackParameters& trackParameters);

};



#endif // #ifndef DDTRACK_CREATOR_ALLEGRO_H
