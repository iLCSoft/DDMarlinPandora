/**
 *  @file   DDMarlinPandora/include/DDGeometryCreatorALLEGRO.h
 * 
 *  @brief  Header file for the geometry creator class.
 * 
 *  $Log: $
 */

#ifndef DDGEOMETRYALLEGRO_CREATOR_H
#define DDGEOMETRYALLEGRO_CREATOR_H 1

#include "Api/PandoraApi.h"

#include "DDRec/DetectorData.h"
#include "DDGeometryCreator.h"

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  DDGeometryCreator class
 */
class DDGeometryCreatorALLEGRO : public DDGeometryCreator
{
public:
    /**
     *  @brief  Constructor
     *
     *  @param  settings the creator settings
     *  @param  pPandora address of the relevant pandora instance
     */
     DDGeometryCreatorALLEGRO(const Settings &settings, const pandora::Pandora *const pPandora);

    /**
     *  @brief  Destructor
     */
     ~DDGeometryCreatorALLEGRO();

    /**
     *  @brief  Create geometry
     */
    pandora::StatusCode CreateGeometry() const;

private:
    /**
     *  @brief  Set mandatory sub detector parameters
     *
     *  @param  subDetectorTypeMap the sub detector type map
     */
    void SetMandatorySubDetectorParameters(SubDetectorTypeMap &subDetectorTypeMap) const;

};

#endif // #ifndef GEOMETRY_CREATOR_H
