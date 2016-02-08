/**
 *  @file   DDMarlinPandora/include/DDTrackCreatorCLIC.h
 * 
 *  @brief  Header file for the ILD implementation of the track creator class.
 * 
 *  $Log: $
 */

#ifndef DDTRACK_CREATOR_CLIC_H
#define DDTRACK_CREATOR_CLIC_H 1

#include "DDTrackCreatorBase.h"

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  DDTrackCreatorCLIC class
 */
class DDTrackCreatorCLIC : public DDTrackCreatorBase
{
public:

    /**
     *  @brief  Constructor
     * 
     *  @param  settings the creator settings
     *  @param  pPandora address of the relevant pandora instance
     */
     DDTrackCreatorCLIC(const Settings &settings, const pandora::Pandora *const pPandora);

    /**
     *  @brief  Destructor
     */
     virtual ~DDTrackCreatorCLIC();

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

};



#endif // #ifndef DDTRACK_CREATOR_CLIC_H
