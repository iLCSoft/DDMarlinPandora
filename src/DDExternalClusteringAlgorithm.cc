/**
 *  @file   DDMarlinPandora/src/DDExternalClusteringAlgorithm.cc
 * 
 *  @brief  Implementation of the external clustering algorithm class.
 * 
 *  $Log: $
 */

#include "EVENT/LCCollection.h"
#include "EVENT/LCEvent.h"
#include "EVENT/Cluster.h"

#include "DDExternalClusteringAlgorithm.h"
#include "DDPandoraPFANewProcessor.h"

#include "Pandora/AlgorithmHeaders.h"

using namespace pandora;

DDExternalClusteringAlgorithm::DDExternalClusteringAlgorithm() :
    m_flagClustersAsPhotons(true)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode DDExternalClusteringAlgorithm::Run()
{
    try
    {
        const CaloHitList *pCaloHitList = NULL;
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pCaloHitList));

        if (pCaloHitList->empty())
            return STATUS_CODE_SUCCESS;

        // Get external photon cluster collection
        const EVENT::LCEvent *const pLCEvent(DDPandoraPFANewProcessor::GetCurrentEvent(&(this->GetPandora())));
        const EVENT::LCCollection *const pExternalClusterCollection = pLCEvent->getCollection(m_externalClusterCollectionName);
        const unsigned int nExternalClusters(pExternalClusterCollection->getNumberOfElements());

        if (0 == nExternalClusters)
            return STATUS_CODE_SUCCESS;

        // Populate pandora parent address to calo hit map
        ParentAddressToCaloHitMap parentAddressToCaloHitMap;

        for (CaloHitList::const_iterator hitIter = pCaloHitList->begin(), hitIterEnd = pCaloHitList->end(); hitIter != hitIterEnd; ++hitIter)
        {
            const pandora::CaloHit *const pCaloHit = *hitIter;
            parentAddressToCaloHitMap.insert(ParentAddressToCaloHitMap::value_type(pCaloHit->GetParentAddress(), pCaloHit));
        }

        // Recreate external clusters within the pandora framework
        for (unsigned int iCluster = 0; iCluster < nExternalClusters; ++iCluster)
        {
            const EVENT::Cluster *const pExternalCluster = dynamic_cast<const EVENT::Cluster*>(pExternalClusterCollection->getElementAt(iCluster));

            if (NULL == pExternalCluster)
                throw EVENT::Exception("Collection type mismatch");

            const CalorimeterHitVec &calorimeterHitVec(pExternalCluster->getCalorimeterHits());

            const pandora::Cluster *pPandoraCluster = NULL;

            for (CalorimeterHitVec::const_iterator iter = calorimeterHitVec.begin(), iterEnd = calorimeterHitVec.end(); iter != iterEnd; ++iter)
            {
                ParentAddressToCaloHitMap::const_iterator pandoraCaloHitIter = parentAddressToCaloHitMap.find(*iter);

                if (parentAddressToCaloHitMap.end() == pandoraCaloHitIter)
                {
                    continue;
                }

                const pandora::CaloHit *const pPandoraCaloHit = pandoraCaloHitIter->second;

                if (NULL == pPandoraCluster)
                {
                    PandoraContentApi::Cluster::Parameters parameters;
                    parameters.m_caloHitList.push_back(pPandoraCaloHit);
                    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Cluster::Create(*this, parameters, pPandoraCluster));

                    if (m_flagClustersAsPhotons)
                    {
                        PandoraContentApi::Cluster::Metadata metadata;
                        metadata.m_particleId = PHOTON;
                        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Cluster::AlterMetadata(*this, pPandoraCluster, metadata));
                    }
                }
                else
                {
                    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AddToCluster(*this, pPandoraCluster, pPandoraCaloHit));
                }
            }
        }
    }
    catch (StatusCodeException &statusCodeException)
    {
        return statusCodeException.GetStatusCode();
    }
    catch (EVENT::Exception &exception)
    {
        std::cout << "DDExternalClusteringAlgorithm failure: " << exception.what() << std::endl;
        return STATUS_CODE_FAILURE;
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode DDExternalClusteringAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle,
        "ExternalClusterCollectionName", m_externalClusterCollectionName));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "FlagClustersAsPhotons", m_flagClustersAsPhotons));

    return STATUS_CODE_SUCCESS;
}
