/**
 *  @file   DDMarlinPandora/include/DDExternalClusteringAlgorithm.h
 * 
 *  @brief  Header file for the external clustering algorithm class.
 * 
 *  $Log: $
 */
#ifndef DDEXTERNALCLUSTERINGALGORITHM_H
#define DDEXTERNALCLUSTERINGALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace pandora { class CaloHit; }

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  DDExternalClusteringAlgorithm class
 */
class DDExternalClusteringAlgorithm : public pandora::Algorithm
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
    DDExternalClusteringAlgorithm();

private:
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    typedef std::map<const void *, const pandora::CaloHit *> ParentAddressToCaloHitMap;

    std::string     m_externalClusterCollectionName = "";   ///< The collection name for the external clusters
    bool            m_flagClustersAsPhotons = true;         ///< Whether to automatically flag new clusters as fixed photons
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *DDExternalClusteringAlgorithm::Factory::CreateAlgorithm() const
{
    return new DDExternalClusteringAlgorithm();
}

#endif // #ifndef DDEXTERNALCLUSTERINGALGORITHM_H
