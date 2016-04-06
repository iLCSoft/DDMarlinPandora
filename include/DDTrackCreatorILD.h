/**
 *  @file   DDMarlinPandora/include/DDTrackCreatorILD.h
 * 
 *  @brief  Header file for the ILD implementation of the track creator class.
 * 
 *  $Log: $
 */

#ifndef DDTRACK_CREATOR_ILD_H
#define DDTRACK_CREATOR_ILD_H 1

#include "DDTrackCreatorBase.h"

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  DDTrackCreatorILD class
 */
class DDTrackCreatorILD : public DDTrackCreatorBase
{
public:

    /**
     *  @brief  Constructor
     * 
     *  @param  settings the creator settings
     *  @param  pPandora address of the relevant pandora instance
     */
     DDTrackCreatorILD(const Settings &settings, const pandora::Pandora *const pPandora);

    /**
     *  @brief  Destructor
     */
     virtual ~DDTrackCreatorILD();

    /**
     *  @brief  Create tracks implementation, insert user code here. Detector specific
     * 
     *  @param  pLCEvent the lcio event
     */
     pandora::StatusCode CreateTracks(EVENT::LCEvent *pLCEvent);
     

    
protected:
    

    //Detector-specific configuration variables
    float                   m_cosTpc;                       ///< Cos(theta) value at end of tpc
    float                   m_tpcInnerR;                    ///< The tpc inner radius
    float                   m_tpcOuterR;                    ///< The tpc outer radius
    unsigned int            m_tpcMaxRow;                    ///< The tpc maximum row number
    float                   m_tpcZmax;                      ///< The tpc maximum z coordinate
    float                   m_tpcMembraneMaxZ;              ///< Tpc membrane max z coordinate

    
    DoubleVector            m_ftdInnerRadii;                ///< List of ftd inner radii
    DoubleVector            m_ftdOuterRadii;                ///< List of ftd outer radii
    DoubleVector            m_ftdZPositions;                ///< List of ftd z positions
    unsigned int            m_nFtdLayers;                   ///< Number of ftd layers
    float                   m_tanLambdaFtd;                 ///< Tan lambda for first ftd layer
    
    float                   m_minEtdZPosition;              ///< Min etd z position
    float                   m_minSetRadius;                 ///< Min set radius

    
    
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
     *  @brief  Get number of hits in TPC of a track
     * 
     *  @param  pTrack the lcio track
     * 
     *  @return number of hits in TPC of a track
     */
    int GetNTpcHits(const EVENT::Track *const pTrack) const;

    /**
     *  @brief  Get number of hits in FTD of a track
     * 
     *  @param  pTrack the lcio track
     * 
     *  @return number of hits in FTDof a track
     */
    int GetNFtdHits(const EVENT::Track *const pTrack) const;

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



#endif // #ifndef DDTRACK_CREATOR_ILD_H
