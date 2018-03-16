/**
 *  @file   DDMarlinPandora/src/DDBFieldPlugin.cc
 * 
 *  @brief  Implementation of the geometry creator class.
 * 
 *  $Log: $
 */

#include "DDBFieldPlugin.h"

#include "Helpers/XmlHelper.h"


DDBFieldPlugin::DDBFieldPlugin(const dd4hep::Detector &detector) :
    m_field(detector.field())
{
    /* nop */
}

//------------------------------------------------------------------------------------------------------------------------------------------

float DDBFieldPlugin::GetBField(const pandora::CartesianVector &positionVector) const
{
    double bfield[3] = {0.};
    m_field.magneticField({positionVector.GetX(), positionVector.GetY(), positionVector.GetZ()}, bfield);
    return std::sqrt(bfield[0]*bfield[0] + bfield[1]*bfield[1] + bfield[2]*bfield[2]);
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode DDBFieldPlugin::Initialize()
{
    /* nop */
    return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode DDBFieldPlugin::ReadSettings(const pandora::TiXmlHandle /*xmlHandle*/)
{
    /* nop */
    return pandora::STATUS_CODE_SUCCESS;
}
