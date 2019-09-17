cmake_minimum_required(VERSION 3.15)

#===================================
#
# SimpleMPM2 Project
#
#===================================
project(SimpleMPM2)

set(VENDORS_DIR ${CMAKE_SOURCE_DIR}/Vendors/)

# GLAD
set(GLAD_INCLUDE_DIR ${VENDORS_DIR}/GLAD/include/)
# GLFW
set(GLFW_INCLUDE_DIR ${VENDORS_DIR}/GLFW/include/)
set(GLFW_LIBRARIES_DIR ${VENDORS_DIR}/GLFW/lib/)
find_library(GLFW_LIBRARIES
    NAMES glfw3
    PATHS ${GLFW_LIBRARIES_DIR}
    )
# glm
set(GLM_INCLUDE_DIR ${VENDORS_DIR}/glm)
# OpenGL headers
set(OPENGL_INCLUDE_DIR
    ${GLAD_INCLUDE_DIR}
    ${GLFW_INCLUDE_DIR}
    ${GLM_INCLUDE_DIR}
    )

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

# Subprojects
# Utilities
add_subdirectory(Utilities/)
# SimulationCore
add_subdirectory(SimulationCore/)
# PostProcessor
add_subdirectory(PostProcessor/)
# TestWithGL (using OpenGL)
add_subdirectory(TestsWithGL/)
