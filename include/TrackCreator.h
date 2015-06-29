/**
 *  @file   MarlinPandora/include/TrackCreator.h
 * 
 *  @brief  Header file for the track creator class.
 * 
 *  $Log: $
 */

#ifndef TRACK_CREATOR_H
#define TRACK_CREATOR_H 1

#include "lcio.h"

#include "IMPL/LCCollectionVec.h"
#include "IMPL/LCFlagImpl.h"

#include "EVENT/LCEvent.h"
#include "EVENT/Track.h"

#include "Api/PandoraApi.h"
#include "Objects/Helix.h"

typedef std::vector<Track *> TrackVector;
typedef std::set<const Track *> TrackList;
typedef std::map<Track *, int> TrackToPidMap;

inline LCCollectionVec *newTrkCol(const std::string &name, LCEvent *evt , bool isSubset)
{
    LCCollectionVec* col = new LCCollectionVec( LCIO::TRACK ) ;

    LCFlagImpl hitFlag(0) ;
    hitFlag.setBit( LCIO::TRBIT_HITS ) ;
    col->setFlag( hitFlag.getFlag()  ) ;
    evt->addCollection( col , name ) ;
    col->setSubset( isSubset ) ;

    return col ;
}

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  TrackCreator class
 */
class TrackCreator
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

        int             m_useOldTrackStateCalculation;          ///< Whether to calculate track states manually, rather than copy stored fitter values

        float           m_d0TrackCut;                           ///< Track d0 cut used to determine whether track can be used to form pfo
        float           m_z0TrackCut;                           ///< Track z0 cut used to determine whether track can be used to form pfo

        int             m_usingNonVertexTracks;                 ///< Whether can form pfos from tracks that don't start at vertex
        int             m_usingUnmatchedNonVertexTracks;        ///< Whether can form pfos from unmatched tracks that don't start at vertex
        int             m_usingUnmatchedVertexTracks;           ///< Whether can form pfos from unmatched tracks that start at vertex
        float           m_unmatchedVertexTrackMaxEnergy;        ///< Maximum energy for unmatched vertex track

        float           m_d0UnmatchedVertexTrackCut;            ///< d0 cut used to determine whether unmatched vertex track can form pfo
        float           m_z0UnmatchedVertexTrackCut;            ///< z0 cut used to determine whether unmatched vertex track can form pfo
        float           m_zCutForNonVertexTracks;               ///< Non vtx track z cut to determine whether track can be used to form pfo

        int             m_reachesECalNTpcHits;                  ///< Minimum number of tpc hits to consider track as reaching ecal
        int             m_reachesECalNFtdHits;                  ///< Minimum number of ftd hits to consider track as reaching ecal
        float           m_reachesECalTpcOuterDistance;          ///< Max distance from track to tpc r max to id whether track reaches ecal
        int             m_reachesECalMinFtdLayer;               ///< Min layer in Ftd for tracks to be considered to have reached decal
        float           m_reachesECalTpcZMaxDistance;           ///< Max distance from track to tpc z max to id whether track reaches ecal
        float           m_reachesECalFtdZMaxDistance;           ///< Max distance from track hit to ftd z position to identify ftd hits
        float           m_curvatureToMomentumFactor;            ///< Constant relating track curvature in b field to momentum

        float           m_minTrackECalDistanceFromIp;           ///< Sanity check on separation between ip and track projected ecal position
        float           m_maxTrackSigmaPOverP;                  ///< Track fraction momentum error cut
        float           m_minMomentumForTrackHitChecks;         ///< Min track momentum required to perform final quality checks on number of hits

        float           m_tpcMembraneMaxZ;                      ///< Tpc membrane max z coordinate
        float           m_maxTpcInnerRDistance;                 ///< Track cut on distance from tpc inner r to id whether track can form pfo
        float           m_minTpcHitFractionOfExpected;          ///< Minimum fraction of TPC hits compared to expected
        int             m_minFtdHitsForTpcHitFraction;          ///< Minimum number of FTD hits to ignore TPC hit fraction
    };

    /**
     *  @brief  Constructor
     * 
     *  @param  settings the creator settings
     *  @param  pPandora address of the relevant pandora instance
     */
     TrackCreator(const Settings &settings, const pandora::Pandora *const pPandora);

    /**
     *  @brief  Destructor
     */
     ~TrackCreator();

    /**
     *  @brief  Create associations between tracks, V0s, kinks, etc
     * 
     *  @param  pLCEvent the lcio event
     */
    pandora::StatusCode CreateTrackAssociations(const EVENT::LCEvent *const pLCEvent);

    /**
     *  @brief  Create tracks, insert user code here
     * 
     *  @param  pLCEvent the lcio event
     */
    pandora::StatusCode CreateTracks(EVENT::LCEvent *pLCEvent);

    /**
     *  @brief  Get the track vector
     * 
     *  @return The track vector
     */
    const TrackVector &GetTrackVector() const;

    /**
     *  @brief  Reset the track creator
     */
    void Reset();

private:
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
    void CopyTrackState(const TrackState *const pTrackState, pandora::InputTrackState &inputTrackState) const;

    /**
     *  @brief  Obtain track states at start and end of track and the momentum at the dca
     * 
     *  @param  pTrack the lcio track
     *  @param  trackParameters the track parameters
     */
    void GetTrackStatesOld(const EVENT::Track *const pTrack, PandoraApi::Track::Parameters &trackParameters) const;

    /**
     *  @brief  Project helix to the surface of the ecal
     * 
     *  @param  pHelix helix fit to be projected to ecal surface
     *  @param  signPz sign w.r.t. increasing z direction
     *  @param  trackParameters the track parameters
     */
    void GetECalProjectionOld(const pandora::Helix *const pHelix, const int signPz, PandoraApi::Track::Parameters &trackParameters) const;

    /**
     *  @brief  Decide whether track reaches the ecal surface
     * 
     *  @param  pTrack the lcio track
     *  @param  trackParameters the track parameters
     */
    void TrackReachesECAL(const EVENT::Track *const pTrack, PandoraApi::Track::Parameters &trackParameters) const;

    /**
     *  @brief  Determine whether a track can be used to form a pfo under the following conditions:
     *          1) if the track proves to be associated with a cluster, OR
     *          2) if the track proves to have no cluster associations
     * 
     *  @param  pTrack the lcio track
     *  @param  trackParameters the track parameters
     */
    void DefineTrackPfoUsage(const EVENT::Track *const pTrack, PandoraApi::Track::Parameters &trackParameters) const;

    /**
     *  @brief  Whether track passes the quality cuts required in order to be used to form a pfo
     * 
     *  @param  pTrack the lcio track
     *  @param  trackParameters the track parameters
     * 
     *  @return boolean
     */
    bool PassesQualityCuts(const EVENT::Track *const pTrack, const PandoraApi::Track::Parameters &trackParameters) const;

    const Settings          m_settings;                     ///< The track creator settings
    const pandora::Pandora *m_pPandora;                     ///< Address of the pandora object to create tracks and track relationships

    const float             m_bField;                       ///< The bfield

    const float             m_tpcInnerR;                    ///< The tpc inner radius
    const float             m_tpcOuterR;                    ///< The tpc outer radius
    const unsigned int      m_tpcMaxRow;                    ///< The tpc maximum row number
    const float             m_tpcZmax;                      ///< The tpc maximum z coordinate
    float                   m_cosTpc;                       ///< Cos(theta) value at end of tpc

    DoubleVector            m_ftdInnerRadii;                ///< List of ftd inner radii
    DoubleVector            m_ftdOuterRadii;                ///< List of ftd outer radii
    DoubleVector            m_ftdZPositions;                ///< List of ftd z positions
    unsigned int            m_nFtdLayers;                   ///< Number of ftd layers
    float                   m_tanLambdaFtd;                 ///< Tan lambda for first ftd layer

    const int               m_eCalBarrelInnerSymmetry;      ///< ECal barrel inner symmetry order
    const float             m_eCalBarrelInnerPhi0;          ///< ECal barrel inner phi 0
    const float             m_eCalBarrelInnerR;             ///< ECal barrel inner radius
    const float             m_eCalEndCapInnerZ;             ///< ECal endcap inner z

    float                   m_minEtdZPosition;              ///< Min etd z position
    float                   m_minSetRadius;                 ///< Min set radius

    TrackVector             m_trackVector;                  ///< The track vector
    TrackList               m_v0TrackList;                  ///< The list of v0 tracks
    TrackList               m_parentTrackList;              ///< The list of parent tracks
    TrackList               m_daughterTrackList;            ///< The list of daughter tracks
    TrackToPidMap           m_trackToPidMap;                ///< The map from track addresses to particle ids, where set by kinks/V0s
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline const TrackVector &TrackCreator::GetTrackVector() const
{
    return m_trackVector;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void TrackCreator::Reset()
{
    m_trackVector.clear();
    m_v0TrackList.clear();
    m_parentTrackList.clear();
    m_daughterTrackList.clear();
    m_trackToPidMap.clear();
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline bool TrackCreator::IsV0(const Track *const pTrack) const
{
    return (m_v0TrackList.end() != m_v0TrackList.find(pTrack));
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline bool TrackCreator::IsParent(const Track *const pTrack) const
{
    return (m_parentTrackList.end() != m_parentTrackList.find(pTrack));
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline bool TrackCreator::IsDaughter(const Track *const pTrack) const
{
    return (m_daughterTrackList.end() != m_daughterTrackList.find(pTrack));
}

#endif // #ifndef TRACK_CREATOR_H
