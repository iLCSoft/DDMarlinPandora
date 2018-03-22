/**
 *  @file   DDMarlinPandora/src/DDBFieldPlugin.cc
 * 
 *  @brief  Implementation of the geometry creator class.
 * 
 *  $Log: $
 */

#include "DDBFieldPlugin.h"

#include "Helpers/XmlHelper.h"

#include "DD4hep/DD4hepUnits.h"

DDBFieldPlugin::DDBFieldPlugin(const dd4hep::Detector &detector) :
    m_field(detector.field())
{
    /* nop */
}

//------------------------------------------------------------------------------------------------------------------------------------------

float DDBFieldPlugin::GetBField(const pandora::CartesianVector &positionVector) const
{
    double bfield[3] = {0.};
    m_field.magneticField({positionVector.GetX()*dd4hep::mm, positionVector.GetY()*dd4hep::mm, positionVector.GetZ()*dd4hep::mm}, bfield);
    return bfield[2]/dd4hep::tesla;
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
