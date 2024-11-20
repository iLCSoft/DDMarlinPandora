########################################################
# cmake file for building DDMarlinPandora
# @author Jan Engels, DESY
CMAKE_MINIMUM_REQUIRED(VERSION 2.8.2 FATAL_ERROR)
########################################################


# project name
PROJECT( DDMarlinPandora )


# project version# continuing from MarlinPandora version 2.0.0 revision 1934
SET( DDMarlinPandora_VERSION_MAJOR 0 )
SET( DDMarlinPandora_VERSION_MINOR 12 )
SET( DDMarlinPandora_VERSION_PATCH 2 )


### CMAKE ###################################################################

# load default settings from ILCSOFT_CMAKE_MODULES
FIND_PACKAGE( ILCUTIL COMPONENTS ILCSOFT_CMAKE_MODULES REQUIRED )

INCLUDE( ilcsoft_default_settings )


### SOURCE ##################################################################

#include directories
INCLUDE_DIRECTORIES( ./include )
INSTALL_DIRECTORY( ./include DESTINATION . FILES_MATCHING PATTERN "*.h" )

AUX_SOURCE_DIRECTORY( ./src/ DDMarlinPandora_SRCS )

### DEPENDENCIES ############################################################

FIND_PACKAGE( Marlin 1.0 REQUIRED )
FIND_PACKAGE( MarlinUtil 01.10 REQUIRED )
FIND_PACKAGE( PandoraSDK 03.00.00 REQUIRED )
FIND_PACKAGE( LCContent 03.00.00 REQUIRED )
FIND_PACKAGE( DD4hep REQUIRED )
FIND_PACKAGE( ROOT REQUIRED )
FIND_PACKAGE( MarlinTrk REQUIRED )

OPTION( USE_SDHCALCONTENT "Use SDHCALContent" OFF )
OPTION( USE_APRILCONTENT "Use APRILContent" OFF )

FOREACH( pkg PandoraSDK Marlin MarlinUtil LCContent DD4hep ROOT MarlinTrk )
    IF( ${pkg}_FOUND )
        INCLUDE_DIRECTORIES( SYSTEM ${${pkg}_INCLUDE_DIRS} )
        LINK_LIBRARIES( ${${pkg}_LIBRARIES} )
        ADD_DEFINITIONS ( ${${pkg}_DEFINITIONS} )
    ENDIF()
ENDFOREACH()

IF( USE_SDHCALCONTENT )
    FIND_PACKAGE( SDHCALContent REQUIRED )
    INCLUDE_DIRECTORIES( SYSTEM ${SDHCALContent_INCLUDE_DIRS} )
    LINK_LIBRARIES( ${SDHCALContent_LIBRARIES} )
    ADD_DEFINITIONS( "-DSDHCALCONTENT" )
ENDIF()

IF( USE_APRILCONTENT )
    SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
    FIND_PACKAGE( APRILContent REQUIRED )
    FIND_PACKAGE(PkgConfig REQUIRED)
    PKG_CHECK_MODULES(MLPACK REQUIRED mlpack)
    INCLUDE_DIRECTORIES(${MLPACK_INCLUDE_DIRS})
    LINK_LIBRARIES(${MLPACK_LIBRARIES})
    INCLUDE_DIRECTORIES( SYSTEM ${APRILContent_INCLUDE_DIRS} )
    LINK_LIBRARIES( ${APRILContent_LIBRARIES} )
    ADD_DEFINITIONS( "-DAPRILCONTENT" )
ENDIF()

IF( PANDORA_MONITORING )
    FIND_PACKAGE( PandoraMonitoring 03.00.00 REQUIRED )
    IF( PandoraMonitoring_FOUND )
        INCLUDE_DIRECTORIES( SYSTEM ${PandoraMonitoring_INCLUDE_DIRS} )
        LINK_LIBRARIES( ${PandoraMonitoring_LIBRARIES} )
        ADD_DEFINITIONS ( ${PandoraMonitoring_DEFINITIONS} )
        ADD_DEFINITIONS( "-DMONITORING" )
    ENDIF()
ENDIF()

### DOCUMENTATION ###########################################################

OPTION( INSTALL_DOC "Set to OFF to skip build/install Documentation" OFF )

IF( INSTALL_DOC AND EXISTS "${PROJECT_SOURCE_DIR}/doc/CMakeLists.txt" )
    ADD_SUBDIRECTORY( ./doc )
ELSE()
    MESSAGE( STATUS "INSTALL_DOC set to OFF" )
    SET( INSTALL_DOC OFF )
ENDIF()


### LIBRARY #################################################################

# add library
ADD_SHARED_LIBRARY( DDMarlinPandora ${DDMarlinPandora_SRCS} )


### INSTALL #################################################################

# install library
INSTALL_SHARED_LIBRARY( DDMarlinPandora DESTINATION lib )

# install header files
INSTALL_DIRECTORY( ./include DESTINATION . )

# display some variables and write them to cache
DISPLAY_STD_VARIABLES()
