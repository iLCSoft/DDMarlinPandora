########################################################
# cmake file for building DDMarlinPandora
# @author Jan Engels, DESY
CMAKE_MINIMUM_REQUIRED(VERSION 2.8.2 FATAL_ERROR)
########################################################


# project name
PROJECT( DDMarlinPandora )


# project version# continuing from MarlinPandora version 2.0.0 revision 1934
SET( DDMarlinPandora_VERSION_MAJOR 00 )
SET( DDMarlinPandora_VERSION_MINOR 01 )
SET( DDMarlinPandora_VERSION_PATCH 00 )


### CMAKE ###################################################################

# load default settings from ILCSOFT_CMAKE_MODULES
FIND_PACKAGE( ILCUTIL COMPONENTS ILCSOFT_CMAKE_MODULES REQUIRED )

INCLUDE( ilcsoft_default_settings )


### SOURCE ##################################################################

#include directories
INCLUDE_DIRECTORIES( ./include )
INSTALL_DIRECTORY( ./include DESTINATION . FILES_MATCHING PATTERN "*.h" )

AUX_SOURCE_DIRECTORY( ./src/ DDMarlinPandora_SRCS )

IF( NOT CMAKE_CXX_FLAGS )
    SET( CMAKE_CXX_FLAGS "-std=c++11" )
    SET( PROCEED_WITHOUT_CXX_FLAGS ON )
ENDIF()

INCLUDE( CheckCXXCompilerFlag )

UNSET( COMPILER_SUPPORTS_CXX_FLAGS CACHE )
CHECK_CXX_COMPILER_FLAG( ${CMAKE_CXX_FLAGS} COMPILER_SUPPORTS_CXX_FLAGS )
IF( NOT COMPILER_SUPPORTS_CXX_FLAGS )
    MESSAGE( STATUS "The compiler ${CMAKE_CXX_COMPILER} does not support cxx flags ${CMAKE_CXX_FLAGS}" )
    IF( NOT PROCEED_WITHOUT_CXX_FLAGS )
        MESSAGE( FATAL_ERROR "CMake abort." )
    ENDIF()
    MESSAGE( STATUS "Proceed without cxx flags ${CMAKE_CXX_FLAGS}" )
    UNSET( CMAKE_CXX_FLAGS )
ENDIF()

if( APPLE )
    SET( CMAKE_CXX_FLAGS "-Wall -pedantic -Wno-long-long -Wno-sign-compare -fno-strict-aliasing ${CMAKE_CXX_FLAGS}" )

else()
    SET( CMAKE_CXX_FLAGS "-Wall -pedantic -Wno-long-long -Wno-sign-compare -fno-strict-aliasing ${CMAKE_CXX_FLAGS}" )
endif()

MESSAGE( STATUS "DDMarlinPandora CMAKE_CXX_FLAGS: ${CMAKE_CXX_FLAGS} " )


### DEPENDENCIES ############################################################

FIND_PACKAGE( Marlin 1.0 REQUIRED )
FIND_PACKAGE( MarlinUtil REQUIRED )
FIND_PACKAGE( PandoraSDK 02.00.00 REQUIRED )
FIND_PACKAGE( LCContent 02.00.00 REQUIRED )
FIND_PACKAGE( DD4hep REQUIRED )
FIND_PACKAGE( ROOT REQUIRED )


FOREACH( pkg Marlin MarlinUtil PandoraSDK LCContent DD4hep ROOT)
    IF( ${pkg}_FOUND )
        INCLUDE_DIRECTORIES( ${${pkg}_INCLUDE_DIRS} )
        LINK_LIBRARIES( ${${pkg}_LIBRARIES} )
        ADD_DEFINITIONS ( ${${pkg}_DEFINITIONS} )
    ENDIF()
ENDFOREACH()

IF( PANDORA_MONITORING )
    FIND_PACKAGE( PandoraMonitoring 02.00.00 REQUIRED )
    IF( PandoraMonitoring_FOUND )
        INCLUDE_DIRECTORIES( ${PandoraMonitoring_INCLUDE_DIRS} )
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
