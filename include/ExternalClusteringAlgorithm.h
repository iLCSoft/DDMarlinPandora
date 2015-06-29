/**
 *  @file   MarlinPandora/include/ExternalClusteringAlgorithm.h
 * 
 *  @brief  Header file for the external clustering algorithm class.
 * 
 *  $Log: $
 */
#ifndef EXTERNAL_CLUSTERING_ALGORITHM_H
#define EXTERNAL_CLUSTERING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace pandora { class CaloHit; }

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  ExternalClusteringAlgorithm class
 */
class ExternalClusteringAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Factory class for instantiating algorithm
     */
    class Factory : public pandora::AlgorithmFactory
    {
    public:
        pandora::Algorithm *CreateAlgorithm() const;
    };

    /**
     *  @brief  Default constructor
     */
    ExternalClusteringAlgorithm();

private:
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    typedef std::map<const void *, const pandora::CaloHit *> ParentAddressToCaloHitMap;

    std::string     m_externalClusterCollectionName;        ///< The collection name for the external clusters
    bool            m_flagClustersAsPhotons;                ///< Whether to automatically flag new clusters as fixed photons
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *ExternalClusteringAlgorithm::Factory::CreateAlgorithm() const
{
    return new ExternalClusteringAlgorithm();
}

#endif // #ifndef EXTERNAL_CLUSTERING_ALGORITHM_H
