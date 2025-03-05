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
    

    //Detector-specific configuration variables
    float                   m_trackerInnerR;                ///< Inner radius of the barrel tracker
    float                   m_trackerOuterR;                ///< Outer radius of the barrel tracker
    float                   m_trackerZmax;                  ///< max extent of the tracker along z
    float                   m_cosTracker;
    
    DoubleVector            m_endcapDiskInnerRadii;                ///< List of endcapDisk inner radii
    DoubleVector            m_endcapDiskOuterRadii;                ///< List of endcapDisk outer radii
    DoubleVector            m_endcapDiskZPositions;                ///< List of endcapDisk z positions
    unsigned int            m_nEndcapDiskLayers;                   ///< Number of endcapDisk layers
    unsigned int            m_barrelTrackerLayers;          ///< Number of barrel tracker layers

    float                   m_tanLambdaEndcapDisk;                 ///< Tan lambda for first endcapDisk layer
    
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
