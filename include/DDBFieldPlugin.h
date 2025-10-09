/**
 *  @file   DDMarlinPandora/include/DDBFieldPlugin.h
 *
 *  @brief  Header file for the BField plugin class.
 *
 *  $Log: $
 */

#ifndef DDBFIELD_PLUGIN_H
#define DDBFIELD_PLUGIN_H 1

#include "Objects/CartesianVector.h"
#include "Plugins/BFieldPlugin.h"

#include "DD4hep/Detector.h"
#include "DD4hep/Fields.h"

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  DDBFieldPlugin class
 */
class DDBFieldPlugin : public pandora::BFieldPlugin {
public:
  DDBFieldPlugin(const dd4hep::Detector& detector);
  float GetBField(const pandora::CartesianVector& positionVector) const;

private:
  pandora::StatusCode Initialize();
  pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

  dd4hep::OverlayedField m_field; ///< The field instance from DD4hep
};

#endif // #ifndef DDBFIELD_PLUGIN_H
