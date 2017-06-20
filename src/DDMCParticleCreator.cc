/**
 *  @file   DDMarlinPandora/src/DDMCParticleCreator.cc
 * 
 *  @brief  Implementation of the mc particle creator class.
 * 
 *  $Log: $
 */

#include "marlin/Global.h"
#include "marlin/Processor.h"

#include "EVENT/LCCollection.h"
#include "EVENT/MCParticle.h"
#include "EVENT/SimCalorimeterHit.h"

#include "UTIL/LCRelationNavigator.h"


#include "DDCaloHitCreator.h"
#include "DDMCParticleCreator.h"
#include "DDTrackCreatorBase.h"

#include <cmath>
#include <limits>

#include "DD4hep/Detector.h"
#include "DD4hep/DD4hepUnits.h"

//forward declarations. See in DDPandoraPFANewProcessor.cc
double getFieldFromCompact();


DDMCParticleCreator::DDMCParticleCreator(const Settings &settings, const pandora::Pandora *const pPandora) :
    m_settings(settings),
    m_pandora(*pPandora),
    m_bField(getFieldFromCompact())
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

DDMCParticleCreator::~DDMCParticleCreator()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode DDMCParticleCreator::CreateMCParticles(const EVENT::LCEvent *const pLCEvent) const
{
    for (StringVector::const_iterator iter = m_settings.m_mcParticleCollections.begin(), iterEnd = m_settings.m_mcParticleCollections.end();
        iter != iterEnd; ++iter)
    {
        try
        {
            const EVENT::LCCollection *pMCParticleCollection = pLCEvent->getCollection(*iter);

            for (int i = 0, iMax = pMCParticleCollection->getNumberOfElements(); i < iMax; ++i)
            {
                try
                {
                    EVENT::MCParticle *pMcParticle = dynamic_cast<MCParticle*>(pMCParticleCollection->getElementAt(i));

                    if (NULL == pMcParticle)
                        throw EVENT::Exception("Collection type mismatch");

                    PandoraApi::MCParticle::Parameters mcParticleParameters;
                    mcParticleParameters.m_energy = pMcParticle->getEnergy();
                    mcParticleParameters.m_particleId = pMcParticle->getPDG();
                    mcParticleParameters.m_mcParticleType = pandora::MC_3D;
                    mcParticleParameters.m_pParentAddress = pMcParticle;
                    mcParticleParameters.m_momentum = pandora::CartesianVector(pMcParticle->getMomentum()[0], pMcParticle->getMomentum()[1],
                        pMcParticle->getMomentum()[2]);
                    mcParticleParameters.m_vertex = pandora::CartesianVector(pMcParticle->getVertex()[0], pMcParticle->getVertex()[1],
                        pMcParticle->getVertex()[2]);
                    mcParticleParameters.m_endpoint = pandora::CartesianVector(pMcParticle->getEndpoint()[0], pMcParticle->getEndpoint()[1],
                        pMcParticle->getEndpoint()[2]);

                    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::MCParticle::Create(m_pandora, mcParticleParameters));

                    // Create parent-daughter relationships
                    for(MCParticleVec::const_iterator itDaughter = pMcParticle->getDaughters().begin(),
                        itDaughterEnd = pMcParticle->getDaughters().end(); itDaughter != itDaughterEnd; ++itDaughter)
                    {
                        PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::SetMCParentDaughterRelationship(m_pandora,
                            pMcParticle, *itDaughter));
                    }
                }
                catch (pandora::StatusCodeException &statusCodeException)
                {
                    streamlog_out(ERROR) << "Failed to extract MCParticle: " << statusCodeException.ToString() << std::endl;
                }
                catch (EVENT::Exception &exception)
                {
                    streamlog_out(WARNING) << "Failed to extract MCParticle: " << exception.what() << std::endl;
                }
            }
        }
        catch (EVENT::Exception &exception)
        {
            streamlog_out(DEBUG5) << "Failed to extract MCParticles collection: " << *iter << ", " << exception.what() << std::endl;
        }
    }

    return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode DDMCParticleCreator::CreateTrackToMCParticleRelationships(const EVENT::LCEvent *const pLCEvent, const TrackVector &trackVector) const
{
    for (StringVector::const_iterator iter = m_settings.m_lcTrackRelationCollections.begin(), iterEnd = m_settings.m_lcTrackRelationCollections.end();
         iter != iterEnd; ++iter)
    {
        try
        {
            const EVENT::LCCollection *pMCRelationCollection = pLCEvent->getCollection(*iter);
            UTIL::LCRelationNavigator navigate(pMCRelationCollection);

            for (TrackVector::const_iterator trackIter = trackVector.begin(), trackIterEnd = trackVector.end();
                trackIter != trackIterEnd; ++trackIter)
            {
                try
                {
                    EVENT::Track *pTrack = *trackIter;
                    const EVENT::LCObjectVec &objectVec = navigate.getRelatedToObjects(*trackIter);

                    // Get reconstructed momentum at dca
                    const pandora::Helix helixFit(pTrack->getPhi(), pTrack->getD0(), pTrack->getZ0(), pTrack->getOmega(), pTrack->getTanLambda(), m_bField);
                    const float recoMomentum(helixFit.GetMomentum().GetMagnitude());

                    // Use momentum magnitude to identify best mc particle
                    MCParticle *pBestMCParticle = NULL;
                    float bestDeltaMomentum(std::numeric_limits<float>::max());

                    for (EVENT::LCObjectVec::const_iterator itRel = objectVec.begin(), itRelEnd = objectVec.end(); itRel != itRelEnd; ++itRel)
                    {
                        EVENT::MCParticle *pMCParticle = NULL;
                        pMCParticle = dynamic_cast<MCParticle *>(*itRel);

                        if (NULL == pMCParticle)
                            continue;

                        const float trueMomentum(pandora::CartesianVector(pMCParticle->getMomentum()[0], pMCParticle->getMomentum()[1],
                            pMCParticle->getMomentum()[2]).GetMagnitude());

                        const float deltaMomentum(std::fabs(recoMomentum - trueMomentum));

                        if (deltaMomentum < bestDeltaMomentum)
                        {
                            pBestMCParticle = pMCParticle;
                            bestDeltaMomentum = deltaMomentum;
                        }
                    }

                    if (NULL == pBestMCParticle)
                        continue;

                    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::SetTrackToMCParticleRelationship(m_pandora, pTrack,
                        pBestMCParticle));
                }
                catch (pandora::StatusCodeException &statusCodeException)
                {
                    streamlog_out(ERROR) << "Failed to extract track to mc particle relationship: " << statusCodeException.ToString() << std::endl;
                }
                catch (EVENT::Exception &exception)
                {
                    streamlog_out(WARNING) << "Failed to extract track to mc particle relationship: " << exception.what() << std::endl;
                }
            }
        }
        catch (EVENT::Exception &exception)
        {
            streamlog_out(DEBUG5) << "Failed to extract track to mc particle relationships collection: " << *iter << ", " << exception.what() << std::endl;
        }
    }

    return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode DDMCParticleCreator::CreateCaloHitToMCParticleRelationships(const EVENT::LCEvent *const pLCEvent, const CalorimeterHitVector &calorimeterHitVector) const
{
    typedef std::map<MCParticle *, float> MCParticleToEnergyWeightMap;
    MCParticleToEnergyWeightMap mcParticleToEnergyWeightMap;

    for (StringVector::const_iterator iter = m_settings.m_lcCaloHitRelationCollections.begin(), iterEnd = m_settings.m_lcCaloHitRelationCollections.end();
         iter != iterEnd; ++iter)
    {
        try
        {
            const EVENT::LCCollection *pMCRelationCollection = pLCEvent->getCollection(*iter);
            UTIL::LCRelationNavigator navigate(pMCRelationCollection);

            for (CalorimeterHitVector::const_iterator caloHitIter = calorimeterHitVector.begin(),
                caloHitIterEnd = calorimeterHitVector.end(); caloHitIter != caloHitIterEnd; ++caloHitIter)
            {
                try
                {
                    mcParticleToEnergyWeightMap.clear();
                    const EVENT::LCObjectVec &objectVec = navigate.getRelatedToObjects(*caloHitIter);

                    for (EVENT::LCObjectVec::const_iterator itRel = objectVec.begin(), itRelEnd = objectVec.end(); itRel != itRelEnd; ++itRel)
                    {
                        EVENT::SimCalorimeterHit *pSimHit = dynamic_cast<SimCalorimeterHit *>(*itRel);

                        if (NULL == pSimHit)
                            continue;

                        for (int iCont = 0, iEnd = pSimHit->getNMCContributions(); iCont < iEnd; ++iCont)
                        {
                            mcParticleToEnergyWeightMap[pSimHit->getParticleCont(iCont)] += pSimHit->getEnergyCont(iCont);
                        }
                    }

                    for (MCParticleToEnergyWeightMap::const_iterator mcParticleIter = mcParticleToEnergyWeightMap.begin(),
                        mcParticleIterEnd = mcParticleToEnergyWeightMap.end(); mcParticleIter != mcParticleIterEnd; ++mcParticleIter)
                    {
                        PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::SetCaloHitToMCParticleRelationship(m_pandora,
                            *caloHitIter, mcParticleIter->first, mcParticleIter->second));
                    }
                }
                catch (pandora::StatusCodeException &statusCodeException)
                {
                    streamlog_out(ERROR) << "Failed to extract calo hit to mc particle relationship: " << statusCodeException.ToString() << std::endl;
                }
                catch (EVENT::Exception &exception)
                {
                    streamlog_out(WARNING) << "Failed to extract calo hit to mc particle relationship: " << exception.what() << std::endl;
                }
            }
        }
        catch (EVENT::Exception &exception)
        {
            streamlog_out(DEBUG5) << "Failed to extract calo hit to mc particle relationships collection: " << *iter << ", " << exception.what() << std::endl;
        }
    }

    return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

DDMCParticleCreator::Settings::Settings():
  m_mcParticleCollections( StringVector() ),
  m_lcCaloHitRelationCollections( StringVector() ),
  m_lcTrackRelationCollections( StringVector() )

{
}
