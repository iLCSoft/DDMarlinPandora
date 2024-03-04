IF(NOT DEFINED mlpack_DIR)
   MESSAGE(STATUS "Warning: it is mandorary to define mlpack_DIR.")
ENDIF()

IF (NOT DEFINED mlpack_INCLUDE_DIRS)
  SET(mlpack_INCLUDE_DIRS ${mlpack_DIR}/include)
ENDIF()

# Check mlpack core header file
FIND_PATH(mlpack_CORE_HPP_DIR NAMES core.hpp HINTS ${mlpack_INCLUDE_DIRS}/mlpack)

# Check mlpack library
#FIND_LIBRARY(mlpack_LIBRARIES NAMES mlpack HINTS ${mlpack_DIR}/*)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(mlpack DEFAULT_MSG mlpack_CORE_HPP_DIR)

IF(NOT mlpack_FOUND)
  MESSAGE(FATAL_ERROR "The mlpack package not found.")
ENDIF()
