cmake_minimum_required(VERSION 3.15)

#===================================
#
# SimpleMPM2 Project
#
#===================================
project(SimpleMPM2)

set(VENDORS_DIR ${CMAKE_SOURCE_DIR}/Vendors/)

#### External dependency
# Eigen
set(EIGEN_INCLUDE_DIR ${VENDORS_DIR}/Eigen/include/)
# GLAD
set(GLAD_INCLUDE_DIR ${VENDORS_DIR}/GLAD/include/)
# GLFW
set(GLFW_INCLUDE_DIR ${VENDORS_DIR}/GLFW/include/)
set(GLFW_LIBRARIES_DIR ${VENDORS_DIR}/GLFW/lib/)
find_library(GLFW_LIBRARIES
    NAMES glfw3
    PATHS ${GLFW_LIBRARIES_DIR}
    )
# OpenGL
set(GLM_INCLUDE_DIR ${VENDORS_DIR}/glm/)
set(OPENGL_INCLUDE_DIR
    ${GLAD_INCLUDE_DIR}
    ${GLFW_INCLUDE_DIR}
    ${GLM_INCLUDE_DIR}
    )
# hdf5
set(HDF5_INCLUDE_DIR ${VENDORS_DIR}/hdf5/include/)
set(HDF5_LIBRARIES_DIR ${VENDORS_DIR}/hdf5/lib/)
#libhdf5 libhdf5_cpp libhdf5_hl libhdf5_hl_cpp libhdf5_tools
find_library(HDF5_LIB
    NAMES libhdf5_D
    PATHS ${HDF5_LIBRARIES_DIR}
    )
find_library(HDF5_CPP_LIB
    NAMES libhdf5_cpp_D
    PATHS ${HDF5_LIBRARIES_DIR}
    )
find_library(HDF5_SZIP
    NAMES libszip_D
    PATHS ${HDF5_LIBRARIES_DIR}
    )
find_library(HDF5_ZLIB
    NAMES libzlib_D
    PATHS ${HDF5_LIBRARIES_DIR}
    )
set(HDF5_LIBRARIES ${HDF5_LIB} ${HDF5_CPP_LIB} ${HDF5_SZIP} ${HDF5_ZLIB})

# Utilities header
set(UTILITIES_INCLUDE_DIR
    ${CMAKE_SOURCE_DIR}/Utilities/
    ${CMAKE_SOURCE_DIR}/Utilities/GraphicUtils
    ${CMAKE_SOURCE_DIR}/Utilities/MemoryUtils
    )

# SimulationCore header
set(SIMULATIONCORE_INCLUDE_DIR
    ${CMAKE_SOURCE_DIR}/SimulationCore/
    ${CMAKE_SOURCE_DIR}/SimulationCore/TriangleMeshRigidBody/
    )

# ConstitutiveModel header
set(CONSTITUTIVEMODEL_INCLUDE_DIR
    ${CMAKE_SOURCE_DIR}/ConstitutiveModel/
    )

# Subprojects
# Utilities
add_subdirectory(Utilities/)
# SimulationCore
add_subdirectory(SimulationCore/)
# ConstitutiveModel
add_subdirectory(ConstitutiveModel/)
# PostProcessor
add_subdirectory(PostProcessor/)
# TestWithGL (using OpenGL)
add_subdirectory(TestsWithGL/)
