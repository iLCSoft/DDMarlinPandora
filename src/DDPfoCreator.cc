/**
 *  @file   DDMarlinPandora/src/DDPfoCreator.cc
 * 
 *  @brief  Implementation of the pfo creator class.
 * 
 *  $Log: $
 */

#include "DDPfoCreator.h"

#include "marlin/Global.h"
#include "marlin/Processor.h"

#include "EVENT/LCCollection.h"

#include "IMPL/ClusterImpl.h"
#include "IMPL/LCCollectionVec.h"
#include "IMPL/LCFlagImpl.h"
#include "IMPL/LCGenericObjectImpl.h"
#include "IMPL/LCRelationImpl.h"
#include "IMPL/ReconstructedParticleImpl.h"
#include "IMPL/VertexImpl.h"

#include "CalorimeterHitType.h"

#include "Api/PandoraApi.h"

#include "Objects/Cluster.h"
#include "Objects/ParticleFlowObject.h"
#include "Objects/Track.h"

#include "Pandora/PdgTable.h"

#include <algorithm>
#include <cmath>

DDPfoCreator::DDPfoCreator(const Settings &settings, const pandora::Pandora *const pPandora) :
    m_settings(settings),
    m_pandora(*pPandora)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

DDPfoCreator::~DDPfoCreator()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode DDPfoCreator::CreateParticleFlowObjects(EVENT::LCEvent *pLCEvent)
{
    const pandora::PfoList *pPandoraPfoList = NULL;
    PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::GetCurrentPfoList(m_pandora, pPandoraPfoList));

    IMPL::LCCollectionVec *pClusterCollection = new IMPL::LCCollectionVec(LCIO::CLUSTER);
    IMPL::LCCollectionVec *pReconstructedParticleCollection = new IMPL::LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);
    IMPL::LCCollectionVec *pStartVertexCollection = new IMPL::LCCollectionVec(LCIO::VERTEX);

    IMPL::LCFlagImpl lcFlagImpl(pClusterCollection->getFlag());
    lcFlagImpl.setBit(LCIO::CLBIT_HITS);
    pClusterCollection->setFlag(lcFlagImpl.getFlag());

    pandora::StringVector subDetectorNames;
    this->InitialiseSubDetectorNames(subDetectorNames);
    pClusterCollection->parameters().setValues("ClusterSubdetectorNames", subDetectorNames);

    // Create lcio "reconstructed particles" from the pandora "particle flow objects"
    for (pandora::PfoList::const_iterator pIter = pPandoraPfoList->begin(), pIterEnd = pPandoraPfoList->end(); pIter != pIterEnd; ++pIter)
    {
        const pandora::ParticleFlowObject *const pPandoraPfo(*pIter);
        IMPL::ReconstructedParticleImpl *const pReconstructedParticle(new ReconstructedParticleImpl());

        const bool hasTrack(!pPandoraPfo->GetTrackList().empty());
        const pandora::ClusterList &clusterList(pPandoraPfo->GetClusterList());

        float clustersTotalEnergy(0.f);
        pandora::CartesianVector referencePoint(0.f, 0.f, 0.f), clustersWeightedPosition(0.f, 0.f, 0.f);
        for (pandora::ClusterList::const_iterator cIter = clusterList.begin(), cIterEnd = clusterList.end(); cIter != cIterEnd; ++cIter)
        {
            const pandora::Cluster *const pPandoraCluster(*cIter);
            pandora::CaloHitList pandoraCaloHitList;
            pPandoraCluster->GetOrderedCaloHitList().FillCaloHitList(pandoraCaloHitList);
            pandoraCaloHitList.insert(pandoraCaloHitList.end(), pPandoraCluster->GetIsolatedCaloHitList().begin(), pPandoraCluster->GetIsolatedCaloHitList().end());

            pandora::FloatVector hitE, hitX, hitY, hitZ;
            IMPL::ClusterImpl *const pLcioCluster(new ClusterImpl());
            this->SetClusterSubDetectorEnergies(subDetectorNames, pLcioCluster, pandoraCaloHitList, hitE, hitX, hitY, hitZ);

            float clusterCorrectEnergy(0.f);
            this->SetClusterEnergyAndError(pPandoraPfo, pPandoraCluster, pLcioCluster, clusterCorrectEnergy);

            pandora::CartesianVector clusterPosition(0.f, 0.f, 0.f);
            const unsigned int nHitsInCluster(pandoraCaloHitList.size());
            this->SetClusterPositionAndError(nHitsInCluster, hitE, hitX, hitY, hitZ, pLcioCluster, clusterPosition);

            if (!hasTrack)
            {
                clustersWeightedPosition += clusterPosition * clusterCorrectEnergy;
                clustersTotalEnergy += clusterCorrectEnergy;
            }

            pClusterCollection->addElement(pLcioCluster);
            pReconstructedParticle->addCluster(pLcioCluster);
        }

        if (!hasTrack)
        {
            if (clustersTotalEnergy < std::numeric_limits<float>::epsilon())
            {
                streamlog_out(WARNING) << "DDPfoCreator::CreateParticleFlowObjects: invalid cluster energy " << clustersTotalEnergy << std::endl;
                throw pandora::StatusCodeException(pandora::STATUS_CODE_FAILURE);
            }
            else
            {
                referencePoint = clustersWeightedPosition * (1.f / clustersTotalEnergy);
            }
        }
        else
        {
            PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, this->CalculateTrackBasedReferencePoint(pPandoraPfo, referencePoint));
        }

        this->SetRecoParticleReferencePoint(referencePoint, pReconstructedParticle);
        this->AddTracksToRecoParticle(pPandoraPfo, pReconstructedParticle);
        this->SetRecoParticlePropertiesFromPFO(pPandoraPfo, pReconstructedParticle);
        pReconstructedParticleCollection->addElement(pReconstructedParticle);

        IMPL::VertexImpl *const pStartVertex(new VertexImpl());
        pStartVertex->setAlgorithmType(m_settings.m_startVertexAlgName.c_str());
        pStartVertex->setPosition(referencePoint.GetX(),referencePoint.GetY(),referencePoint.GetZ());
        pStartVertex->setAssociatedParticle(pReconstructedParticle);
        pStartVertexCollection->addElement(pStartVertex);
    }

    pLCEvent->addCollection(pClusterCollection, m_settings.m_clusterCollectionName.c_str());
    pLCEvent->addCollection(pReconstructedParticleCollection, m_settings.m_pfoCollectionName.c_str());
    pLCEvent->addCollection(pStartVertexCollection, m_settings.m_startVertexCollectionName.c_str());

    return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DDPfoCreator::InitialiseSubDetectorNames(pandora::StringVector &subDetectorNames) const
{
    subDetectorNames.push_back("ecal");
    subDetectorNames.push_back("hcal");
    subDetectorNames.push_back("yoke");
    subDetectorNames.push_back("lcal");
    subDetectorNames.push_back("lhcal");
    subDetectorNames.push_back("bcal");
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DDPfoCreator::SetClusterSubDetectorEnergies(const pandora::StringVector &subDetectorNames, IMPL::ClusterImpl *const pLcioCluster,
    const pandora::CaloHitList &pandoraCaloHitList, pandora::FloatVector &hitE, pandora::FloatVector &hitX, pandora::FloatVector &hitY,
    pandora::FloatVector &hitZ) const
{
    for (pandora::CaloHitList::const_iterator hIter = pandoraCaloHitList.begin(), hIterEnd = pandoraCaloHitList.end(); hIter != hIterEnd; ++hIter)
    {
        const pandora::CaloHit *const pPandoraCaloHit(*hIter);
        EVENT::CalorimeterHit *const pCalorimeterHit = (EVENT::CalorimeterHit*)(pPandoraCaloHit->GetParentAddress());
        pLcioCluster->addHit(pCalorimeterHit, 1.f);

        const float caloHitEnergy(pCalorimeterHit->getEnergy());
        hitE.push_back(caloHitEnergy);
        hitX.push_back(pCalorimeterHit->getPosition()[0]);
        hitY.push_back(pCalorimeterHit->getPosition()[1]);
        hitZ.push_back(pCalorimeterHit->getPosition()[2]);

        std::vector<float> &subDetectorEnergies = pLcioCluster->subdetectorEnergies();
        subDetectorEnergies.resize(subDetectorNames.size());

        switch (CHT(pCalorimeterHit->getType()).caloID())
        {
            case CHT::ecal:  subDetectorEnergies[ECAL_INDEX ] += caloHitEnergy; break;
            case CHT::hcal:  subDetectorEnergies[HCAL_INDEX ] += caloHitEnergy; break;
            case CHT::yoke:  subDetectorEnergies[YOKE_INDEX ] += caloHitEnergy; break;
            case CHT::lcal:  subDetectorEnergies[LCAL_INDEX ] += caloHitEnergy; break;
            case CHT::lhcal: subDetectorEnergies[LHCAL_INDEX] += caloHitEnergy; break;
            case CHT::bcal:  subDetectorEnergies[BCAL_INDEX ] += caloHitEnergy; break;
            default: streamlog_out(WARNING) << "DDPfoCreator::SetClusterSubDetectorEnergies: no subdetector found for hit with type: " << pCalorimeterHit->getType() << std::endl;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DDPfoCreator::SetClusterEnergyAndError(const pandora::ParticleFlowObject *const pPandoraPfo, const pandora::Cluster *const pPandoraCluster, 
    IMPL::ClusterImpl *const pLcioCluster, float &clusterCorrectEnergy) const
{
    const bool isEmShower((pandora::PHOTON == pPandoraPfo->GetParticleId()) || (pandora::E_MINUS == std::abs(pPandoraPfo->GetParticleId())));
    clusterCorrectEnergy = (isEmShower ? pPandoraCluster->GetCorrectedElectromagneticEnergy(m_pandora) : pPandoraCluster->GetCorrectedHadronicEnergy(m_pandora));

    if (clusterCorrectEnergy < std::numeric_limits<float>::epsilon())
        throw pandora::StatusCodeException(pandora::STATUS_CODE_FAILURE);

    const float stochasticTerm(isEmShower ? m_settings.m_emStochasticTerm : m_settings.m_hadStochasticTerm); 
    const float constantTerm(isEmShower ? m_settings.m_emConstantTerm : m_settings.m_hadConstantTerm);
    const float energyError(std::sqrt(stochasticTerm * stochasticTerm / clusterCorrectEnergy + constantTerm * constantTerm) * clusterCorrectEnergy);

    pLcioCluster->setEnergy(clusterCorrectEnergy);
    pLcioCluster->setEnergyError(energyError);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DDPfoCreator::SetClusterPositionAndError(const unsigned int nHitsInCluster, pandora::FloatVector &hitE, pandora::FloatVector &hitX, 
    pandora::FloatVector &hitY, pandora::FloatVector &hitZ, IMPL::ClusterImpl *const pLcioCluster, pandora::CartesianVector &clusterPositionVec) const
{
    ClusterShapes *const pClusterShapes(new ClusterShapes(nHitsInCluster, hitE.data(), hitX.data(), hitY.data(), hitZ.data()));

    try
    {
        pLcioCluster->setIPhi(std::atan2(pClusterShapes->getEigenVecInertia()[1], pClusterShapes->getEigenVecInertia()[0]));
        pLcioCluster->setITheta(std::acos(pClusterShapes->getEigenVecInertia()[2]));
        pLcioCluster->setPosition(pClusterShapes->getCentreOfGravity());
        //ATTN these two lines below would only compile with ilcsoft HEAD V2015-10-13 and above
        //pLcioCluster->setPositionError(pClusterShapes->getCenterOfGravityErrors());
        //pLcioCluster->setDirectionError(pClusterShapes->getEigenVecInertiaErrors());
        clusterPositionVec.SetValues(pClusterShapes->getCentreOfGravity()[0], pClusterShapes->getCentreOfGravity()[1], pClusterShapes->getCentreOfGravity()[2]);
    }
    catch (...)
    {
        streamlog_out(WARNING) << "DDPfoCreator::SetClusterPositionAndError: unidentified exception caught." << std::endl;
    }

    delete pClusterShapes;
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode DDPfoCreator::CalculateTrackBasedReferencePoint(const pandora::ParticleFlowObject *const pPandoraPfo, pandora::CartesianVector &referencePoint) const
{
    const pandora::TrackList &trackList(pPandoraPfo->GetTrackList());

    float totalTrackMomentumAtDca(0.f), totalTrackMomentumAtStart(0.f);
    pandora::CartesianVector referencePointAtDCAWeighted(0.f, 0.f, 0.f), referencePointAtStartWeighted(0.f, 0.f, 0.f);

    bool hasSiblings(false);
    for (pandora::TrackList::const_iterator tIter = trackList.begin(), tIterEnd = trackList.end(); tIter != tIterEnd; ++tIter)
    {
        const pandora::Track *const pPandoraTrack(*tIter);

        if (!this->IsValidParentTrack(pPandoraTrack, trackList))
            continue;

        if (this->HasValidSiblingTrack(pPandoraTrack, trackList))
        {
            // Presence of sibling tracks typically represents a conversion
            const pandora::CartesianVector &trackStartPoint((pPandoraTrack->GetTrackStateAtStart()).GetPosition());
            const float trackStartMomentum(((pPandoraTrack->GetTrackStateAtStart()).GetMomentum()).GetMagnitude());
            referencePointAtStartWeighted += trackStartPoint * trackStartMomentum;
            totalTrackMomentumAtStart += trackStartMomentum;
            hasSiblings = true;
        }
        else
        {
            const EVENT::Track *const pLcioTrack = (EVENT::Track*)(pPandoraTrack->GetParentAddress());
            const float z0(pPandoraTrack->GetZ0());
            pandora::CartesianVector intersectionPoint(0.f, 0.f, 0.f);

            intersectionPoint.SetValues(pLcioTrack->getD0() * std::cos(pLcioTrack->getPhi()), pLcioTrack->getD0() * std::sin(pLcioTrack->getPhi()), z0);
            const float trackMomentumAtDca((pPandoraTrack->GetMomentumAtDca()).GetMagnitude());
            referencePointAtDCAWeighted += intersectionPoint * trackMomentumAtDca;
            totalTrackMomentumAtDca += trackMomentumAtDca;
        }
    }

    if (hasSiblings)
    {
        if (totalTrackMomentumAtStart < std::numeric_limits<float>::epsilon())
        {
            streamlog_out(WARNING) << "DDPfoCreator::CalculateTrackBasedReferencePoint: invalid track momentum " << totalTrackMomentumAtStart << std::endl;
            throw pandora::StatusCodeException(pandora::STATUS_CODE_FAILURE);
        }
        else
        {
            referencePoint = referencePointAtStartWeighted * (1.f / totalTrackMomentumAtStart);
        }
    }
    else
    {
        if (totalTrackMomentumAtDca < std::numeric_limits<float>::epsilon())
        {
            streamlog_out(WARNING) << "DDPfoCreator::CalculateTrackBasedReferencePoint: invalid track momentum " << totalTrackMomentumAtDca << std::endl;
            throw pandora::StatusCodeException(pandora::STATUS_CODE_FAILURE);
        }
        else
        {
            referencePoint = referencePointAtDCAWeighted * (1.f / totalTrackMomentumAtDca);
        }
    }

    return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool DDPfoCreator::IsValidParentTrack(const pandora::Track *const pPandoraTrack, const pandora::TrackList &allTrackList) const
{
    const pandora::TrackList &parentTrackList(pPandoraTrack->GetParentList());

    for (pandora::TrackList::const_iterator iter = parentTrackList.begin(), iterEnd = parentTrackList.end(); iter != iterEnd; ++iter)
    {
        if (allTrackList.end() != std::find(allTrackList.begin(), allTrackList.end(), *iter))
            continue;

        // ATTN This track must have a parent not in the all track list; still use it if it is the closest to the ip
        streamlog_out(WARNING) << "DDPfoCreator::IsValidParentTrack: mismatch in track relationship information, use information as available " << std::endl;

        if (this->IsClosestTrackToIP(pPandoraTrack, allTrackList))
            return true;

        return false;
    }

    // Ideal case: All parents are associated to same pfo
    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool DDPfoCreator::HasValidSiblingTrack(const pandora::Track *const pPandoraTrack, const pandora::TrackList &allTrackList) const
{
    const pandora::TrackList &siblingTrackList(pPandoraTrack->GetSiblingList());

    for (pandora::TrackList::const_iterator iter = siblingTrackList.begin(), iterEnd = siblingTrackList.end(); iter != iterEnd; ++iter)
    {
        if (allTrackList.end() != std::find(allTrackList.begin(), allTrackList.end(), *iter))
            continue;

        // ATTN This track must have a sibling not in the all track list; still use it if it has a second sibling that is in the list
        streamlog_out(WARNING) << "DDPfoCreator::HasValidSiblingTrack: mismatch in track relationship information, use information as available " << std::endl;

        if (this->AreAnyOtherSiblingsInList(pPandoraTrack, allTrackList))
            return true;

        return false;
    }

    // Ideal case: All siblings associated to same pfo
    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool DDPfoCreator::IsClosestTrackToIP(const pandora::Track *const pPandoraTrack, const pandora::TrackList &allTrackList) const 
{
    const pandora::Track *pClosestTrack(NULL);
    float closestTrackDisplacement(std::numeric_limits<float>::max()); 

    for (pandora::TrackList::const_iterator iter = allTrackList.begin(), iterEnd = allTrackList.end(); iter != iterEnd; ++iter)
    {
        const pandora::Track *const pTrack(*iter);
        const float trialTrackDisplacement(pTrack->GetTrackStateAtStart().GetPosition().GetMagnitude());

        if (trialTrackDisplacement < closestTrackDisplacement)
        {
            closestTrackDisplacement = trialTrackDisplacement;
            pClosestTrack = pTrack;
        }
    }

    return (pPandoraTrack == pClosestTrack);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool DDPfoCreator::AreAnyOtherSiblingsInList(const pandora::Track *const pPandoraTrack, const pandora::TrackList &allTrackList) const
{
    const pandora::TrackList &siblingTrackList(pPandoraTrack->GetSiblingList());

    for (pandora::TrackList::const_iterator iter = siblingTrackList.begin(), iterEnd = siblingTrackList.end(); iter != iterEnd; ++iter)
    {
        if (allTrackList.end() != std::find(allTrackList.begin(), allTrackList.end(), *iter))
            return true;
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DDPfoCreator::SetRecoParticleReferencePoint(const pandora::CartesianVector &referencePoint, IMPL::ReconstructedParticleImpl *const pReconstructedParticle) const
{
    const float referencePointArray[3] = {referencePoint.GetX(), referencePoint.GetY(), referencePoint.GetZ()};
    pReconstructedParticle->setReferencePoint(referencePointArray);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DDPfoCreator::AddTracksToRecoParticle(const pandora::ParticleFlowObject *const pPandoraPfo, IMPL::ReconstructedParticleImpl *const pReconstructedParticle) const
{
    const pandora::TrackList &trackList(pPandoraPfo->GetTrackList());

    for (pandora::TrackList::const_iterator tIter = trackList.begin(), tIterEnd = trackList.end(); tIter != tIterEnd; ++tIter)
    {
        const pandora::Track *const pTrack(*tIter);
        pReconstructedParticle->addTrack((EVENT::Track*)(pTrack->GetParentAddress()));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DDPfoCreator::SetRecoParticlePropertiesFromPFO(const pandora::ParticleFlowObject *const pPandoraPfo, IMPL::ReconstructedParticleImpl *const pReconstructedParticle) const
{
    const float momentum[3] = {pPandoraPfo->GetMomentum().GetX(), pPandoraPfo->GetMomentum().GetY(), pPandoraPfo->GetMomentum().GetZ()};
    pReconstructedParticle->setMomentum(momentum);
    pReconstructedParticle->setEnergy(pPandoraPfo->GetEnergy());
    pReconstructedParticle->setMass(pPandoraPfo->GetMass());
    pReconstructedParticle->setCharge(pPandoraPfo->GetCharge());
    pReconstructedParticle->setType(pPandoraPfo->GetParticleId());
}

//------------------------------------------------------------------------------------------------------------------------------------------

DDPfoCreator::Settings::Settings():
    m_clusterCollectionName(""),
    m_pfoCollectionName(""),
    m_startVertexCollectionName(""),
    m_startVertexAlgName(""),
    m_emStochasticTerm(0.17f),
    m_hadStochasticTerm(0.6f),
    m_emConstantTerm(0.01f),
    m_hadConstantTerm(0.03f)
{
}
