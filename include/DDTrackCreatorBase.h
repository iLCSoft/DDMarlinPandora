/**
 *  @file   DDMarlinPandora/include/DDTrackCreatorBase.h
 * 
 *  @brief  Header file for the track creator base class.
 * 
 *  $Log: $
 */

#ifndef DDTRACK_CREATOR_BASE_H
#define DDTRACK_CREATOR_BASE_H 1

#include "lcio.h"

#include "IMPL/LCCollectionVec.h"
#include "IMPL/LCFlagImpl.h"

#include "EVENT/LCEvent.h"
#include "EVENT/Track.h"

#include "Api/PandoraApi.h"
#include "Objects/Helix.h"

#include <UTIL/ILDConf.h>

#include <memory>

typedef std::vector<EVENT::Track *> TrackVector;
typedef std::set<const EVENT::Track *> TrackList;
typedef std::map<EVENT::Track *, int> TrackToPidMap;

namespace lc_content{
  class LCTrackParameters;
  class LCTrackFactory;
}
namespace MarlinTrk{
  class IMarlinTrkSystem;
}

inline IMPL::LCCollectionVec *newTrkCol(const std::string &name, EVENT::LCEvent *evt , bool isSubset)
{
    IMPL::LCCollectionVec* col = new IMPL::LCCollectionVec( EVENT::LCIO::TRACK ) ;

    IMPL::LCFlagImpl hitFlag(0) ;
    hitFlag.setBit( EVENT::LCIO::TRBIT_HITS ) ;
    col->setFlag( hitFlag.getFlag()  ) ;
    evt->addCollection( col , name ) ;
    col->setSubset( isSubset ) ;

    return col ;
}

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  DDTrackCreatorBase class
 */
class DDTrackCreatorBase
{
public:
    typedef std::vector<double> DoubleVector;
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

        StringVector    m_trackCollections;                     ///< The reconstructed track collections
        StringVector    m_kinkVertexCollections;                ///< The kink vertex collections
        StringVector    m_prongVertexCollections;               ///< The prong vertex collections
        StringVector    m_splitVertexCollections;               ///< The split vertex collections
        StringVector    m_v0VertexCollections;                  ///< The v0 vertex collections

        StringVector    m_prongSplitVertexCollections;          ///< Concatenated list of prong and split vertex collections
        int             m_shouldFormTrackRelationships;         ///< Whether to form pandora track relationships using v0 and kink info

        int             m_minTrackHits;                         ///< Track quality cut: the minimum number of track hits
        int             m_minFtdTrackHits;                      ///< Track quality cut: the minimum number of FTD track hits for FTD only tracks
        int             m_maxTrackHits;                         ///< Track quality cut: the maximum number of track hits

        float           m_d0TrackCut;                           ///< Track d0 cut used to determine whether track can be used to form pfo
        float           m_z0TrackCut;                           ///< Track z0 cut used to determine whether track can be used to form pfo

        int             m_usingNonVertexTracks;                 ///< Whether can form pfos from tracks that don't start at vertex
        int             m_usingUnmatchedNonVertexTracks;        ///< Whether can form pfos from unmatched tracks that don't start at vertex
        int             m_usingUnmatchedVertexTracks;           ///< Whether can form pfos from unmatched tracks that start at vertex
        float           m_unmatchedVertexTrackMaxEnergy;        ///< Maximum energy for unmatched vertex track

        float           m_d0UnmatchedVertexTrackCut;            ///< d0 cut used to determine whether unmatched vertex track can form pfo
        float           m_z0UnmatchedVertexTrackCut;            ///< z0 cut used to determine whether unmatched vertex track can form pfo
        float           m_zCutForNonVertexTracks;               ///< Non vtx track z cut to determine whether track can be used to form pfo

        int             m_reachesECalNBarrelTrackerHits;                  ///< Minimum number of barrel tracker hits to consider track as reaching ecal
        int             m_reachesECalNFtdHits;                  ///< Minimum number of ftd hits to consider track as reaching ecal
        float           m_reachesECalBarrelTrackerOuterDistance;          ///< Max distance from track to barrel tracker r max to id whether track reaches ecal
        int             m_reachesECalMinFtdLayer;               ///< Min layer in Ftd for tracks to be considered to have reached decal
        float           m_reachesECalBarrelTrackerZMaxDistance;           ///< Max distance from track to barrel tracker z max to id whether track reaches ecal
        float           m_reachesECalFtdZMaxDistance;           ///< Max distance from track hit to ftd z position to identify ftd hits
        float           m_curvatureToMomentumFactor;            ///< Constant relating track curvature in b field to momentum

        float           m_minTrackECalDistanceFromIp;           ///< Sanity check on separation between ip and track projected ecal position
        float           m_maxTrackSigmaPOverP;                  ///< Track fraction momentum error cut
        float           m_minMomentumForTrackHitChecks;         ///< Min track momentum required to perform final quality checks on number of hits

        float           m_maxBarrelTrackerInnerRDistance;                 ///< Track cut on distance from barrel tracker inner r to id whether track can form pfo
        float           m_minBarrelTrackerHitFractionOfExpected;          ///< Minimum fraction of TPC hits compared to expected
        int             m_minFtdHitsForBarrelTrackerHitFraction;          ///< Minimum number of FTD hits to ignore TPC hit fraction
        float           m_trackStateTolerance; ///< distance below tracker ecal radius the second trackstate in the ecal endcap is still passed to pandora
        std::string     m_trackingSystemName;  ///< name of the tracking system used for getting new track states


        ///Nikiforos: Moved from main class
        
        float             m_bField;                       ///< The bfield
        int               m_eCalBarrelInnerSymmetry;      ///< ECal barrel inner symmetry order
        float             m_eCalBarrelInnerPhi0;          ///< ECal barrel inner phi 0
        float             m_eCalBarrelInnerR;             ///< ECal barrel inner radius
        float             m_eCalEndCapInnerZ;             ///< ECal endcap inner z

        //Tracking Detector names not needed anymore, accessed by det type flag

        
    };

    /**
     *  @brief  Constructor
     * 
     *  @param  settings the creator settings
     *  @param  pPandora address of the relevant pandora instance
     */
     DDTrackCreatorBase(const Settings &settings, const pandora::Pandora *const pPandora);

    /**
     *  @brief  Destructor
     */
     virtual ~DDTrackCreatorBase();

    /**
     *  @brief  Create associations between tracks, V0s, kinks, etc
     * 
     *  @param  pLCEvent the lcio event
     */
    pandora::StatusCode CreateTrackAssociations(const EVENT::LCEvent *const pLCEvent);

    /**
     *  @brief  Create tracks, insert user code here. Implement accordin to detector model
     * 
     *  @param  pLCEvent the lcio event
     */
    virtual pandora::StatusCode CreateTracks(EVENT::LCEvent *pLCEvent) = 0 ;

    /**
     *  @brief  Get the track vector
     * 
     *  @return The track vector
     */
    const TrackVector &GetTrackVector() const;


    /**
     *  @brief  Calculate possible second track state at the ECal Endcap
     *
     *  @param track lcio track
     *  @param trackParameters pandora LCTrackParameters
     */
    virtual void GetTrackStatesAtCalo( EVENT::Track *track, lc_content::LCTrackParameters& trackParameters);


    /**
     *  @brief  Reset the track creator
     */
    void Reset();
    
    
    

protected:
    
    
    

    const Settings          m_settings;                     ///< The track creator settings
    const pandora::Pandora &m_pandora;                     ///< Reference to the pandora object to create tracks and track relationships


    TrackVector             m_trackVector;                  ///< The track vector
    TrackList               m_v0TrackList;                  ///< The list of v0 tracks
    TrackList               m_parentTrackList;              ///< The list of parent tracks
    TrackList               m_daughterTrackList;            ///< The list of daughter tracks
    TrackToPidMap           m_trackToPidMap;                ///< The map from track addresses to particle ids, where set by kinks/V0s
    float                   m_minimalTrackStateRadiusSquared;         ///< minimal track state radius, derived value
    std::shared_ptr<MarlinTrk::IMarlinTrkSystem> m_trackingSystem={}; ///< Tracking system used for track states
    std::shared_ptr<UTIL::BitField64> m_encoder={};                   ///< cell ID encoder
    std::shared_ptr<lc_content::LCTrackFactory> m_lcTrackFactory={};  ///< LCTrackFactor for creating LCTracks

    
    ///Nikiforos: Need to implement following abstract functions according to detector model
    
     /**
     *  @brief  Whether track passes the quality cuts required in order to be used to form a pfo
     * 
     *  @param  pTrack the lcio track
     *  @param  trackParameters the track parameters
     * 
     *  @return boolean
     */
    virtual bool PassesQualityCuts(const EVENT::Track *const pTrack, const PandoraApi::Track::Parameters &trackParameters) const = 0;
    
    /**
     *  @brief  Decide whether track reaches the ecal surface
     * 
     *  @param  pTrack the lcio track
     *  @param  trackParameters the track parameters
     */
    virtual void TrackReachesECAL(const EVENT::Track *const pTrack, PandoraApi::Track::Parameters &trackParameters) const = 0 ;

    /**
     *  @brief  Determine whether a track can be used to form a pfo under the following conditions:
     *          1) if the track proves to be associated with a cluster, OR
     *          2) if the track proves to have no cluster associations
     * 
     *  @param  pTrack the lcio track
     *  @param  trackParameters the track parameters
     */
    virtual void DefineTrackPfoUsage(const EVENT::Track *const pTrack, PandoraApi::Track::Parameters &trackParameters) const = 0;
    
    /**
     *  @brief  Extract kink information from specified lcio collections
     * 
     *  @param  pLCEvent the lcio event
     */
    pandora::StatusCode ExtractKinks(const EVENT::LCEvent *const pLCEvent);

    /**
     *  @brief  Extract prong and split information from specified lcio collections
     * 
     *  @param  pLCEvent the lcio event
     */
    pandora::StatusCode ExtractProngsAndSplits(const EVENT::LCEvent *const pLCEvent);

    /**
     *  @brief  Extract v0 information from specified lcio collections
     * 
     *  @param  pLCEvent the lcio event
     */
    pandora::StatusCode ExtractV0s(const EVENT::LCEvent *const pLCEvent);

    /**
     *  @brief  Whether the track vertex conflicts with previously provided relationship information
     * 
     *  @param  trackVec the vector of tracks associated with the vertex
     */
    bool IsConflictingRelationship(const EVENT::TrackVec &trackVec) const;

    /**
     *  @brief  Whether a track is a v0 track
     * 
     *  @param  pTrack the lcio track
     * 
     *  @return boolean
     */
    bool IsV0(const EVENT::Track *const pTrack) const;

    /**
     *  @brief  Whether a track is a parent track
     * 
     *  @param  pTrack the lcio track
     * 
     *  @return boolean
     */
    bool IsParent(const EVENT::Track *const pTrack) const;

    /**
     *  @brief  Whether a track is a daughter track
     * 
     *  @param  pTrack the lcio track
     * 
     *  @return boolean
     */
    bool IsDaughter(const EVENT::Track *const pTrack) const;

    /**
     *  @brief  Copy track states stored in lcio tracks to pandora track parameters
     * 
     *  @param  pTrack the lcio track
     *  @param  trackParameters the track parameters
     */
    void GetTrackStates(const EVENT::Track *const pTrack, PandoraApi::Track::Parameters &trackParameters) const;

    /**
     *  @brief  Copy track state from lcio track state instance to pandora input track state
     * 
     *  @param  pTrackState the lcio track state instance
     *  @param  inputTrackState the pandora input track state
     */
    void CopyTrackState(const EVENT::TrackState *const pTrackState, pandora::InputTrackState &inputTrackState) const;

    /**
     *  @brief  Obtain track time when it reaches ECAL
     * 
     *  @param  pTrack the lcio track
     */
    float CalculateTrackTimeAtCalorimeter(const EVENT::Track *const pTrack) const;
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline const TrackVector &DDTrackCreatorBase::GetTrackVector() const
{
    return m_trackVector;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void DDTrackCreatorBase::Reset()
{
    m_trackVector.clear();
    m_v0TrackList.clear();
    m_parentTrackList.clear();
    m_daughterTrackList.clear();
    m_trackToPidMap.clear();
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline bool DDTrackCreatorBase::IsV0(const EVENT::Track *const pTrack) const
{
    return (m_v0TrackList.end() != m_v0TrackList.find(pTrack));
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline bool DDTrackCreatorBase::IsParent(const EVENT::Track *const pTrack) const
{
    return (m_parentTrackList.end() != m_parentTrackList.find(pTrack));
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline bool DDTrackCreatorBase::IsDaughter(const EVENT::Track *const pTrack) const
{
    return (m_daughterTrackList.end() != m_daughterTrackList.find(pTrack));
}

#endif // #ifndef DDTRACK_CREATOR_BASE_H
