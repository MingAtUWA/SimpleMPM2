#===================================
#
# SimulationCore Project
#
#===================================
add_library(
    SimulationCore STATIC
    SimulationCore_pcp.h SimulationCore_pcp.cpp
    BC.h
    # Base Class
    Model.h
    Step.h Step.cpp
    TimeHistory.h
    ResultFile.h
    # ProgressBar
    TimeHistory_ConsoleProgressBar.h TimeHistory_ConsoleProgressBar.cpp
    # Plain binary result file
    ResultFile_PlainBin.h ResultFile_PlainBin.cpp
    ResultFile_XML.h ResultFile_XML.cpp
    # Mechanics, material point and rigid-body smooth contact
    Model_S2D_ME_s_RigidBody.h Model_S2D_ME_s_RigidBody.cpp
    Step_S2D_ME_s_RigidBody.h Step_S2D_ME_s_RigidBody.cpp
    TimeHistory_S2D_ME_s_RigidBody.h
    #
    Step_S2D_ME_s_RigidBody_Fric.h Step_S2D_ME_s_RigidBody_Fric.cpp
    TimeHistory_S2D_ME_s_RigidBody_Fric.h
    )

add_subdirectory(TriangleMeshRigidBody/)

target_include_directories(
    SimulationCore
    PRIVATE
        ${SIMULATIONCORE_INCLUDE_DIR}
        ${UTILITIES_INCLUDE_DIR}
        ${OPENGL_INCLUDE_DIR}
    )

target_link_libraries(
    SimulationCore
    PRIVATE
        Utilities
    )

if (MSVC)
    target_compile_definitions(SimulationCore PRIVATE _CRT_SECURE_NO_WARNINGS)
    target_compile_options(SimulationCore PRIVATE -YcSimulationCore_pcp.h)
endif()
