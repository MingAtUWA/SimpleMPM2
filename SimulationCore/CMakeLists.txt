#===================================
#
# SimulationCore Project
#
#===================================
add_library(
    SimulationCore STATIC
    SimulationCore_pcp.h SimulationCore_pcp.cpp
    BC.h
    # Utilities
    TriangleMeshToParticles.h TriangleMeshToParticles.cpp
    # Base Class
    Model.h
    Step.h Step.cpp
    ModelDataOutput.h ModelDataOutput.cpp
    TimeHistoryOutput.h TimeHistoryOutput.cpp
    ResultFile.h
    # Plain binary result file
    ResultFile_PlainBin.h ResultFile_PlainBin_DataStruct.h
    # XML result file
    ResultFile_XML.h
    # ProgressBar
    TimeHistoryOutput_ConsoleProgressBar.h TimeHistoryOutput_ConsoleProgressBar.cpp
    # Purely mechanics, material point and rigid-body smooth contact
    Model_S2D_ME_s_RigidBody.h Model_S2D_ME_s_RigidBody.cpp
    Step_S2D_ME_s_RigidBody.h Step_S2D_ME_s_RigidBody.cpp
    ModelDataOutput_S2D_ME_s_RigidBody.h ModelDataOutput_S2D_ME_s_RigidBody.cpp
    TimeHistoryOutput_S2D_ME_s_RigidBody.h  TimeHistoryOutput_S2D_ME_s_RigidBody.cpp
    # Purely mechanics, material ponit and rigid-body frictional contact
    ContactState.h # keep tract of contact state
    Model_S2D_ME_s_RigidBody_Fric.h Model_S2D_ME_s_RigidBody_Fric.cpp
    Step_S2D_ME_s_RigidBody_Fric.h Step_S2D_ME_s_RigidBody_Fric.cpp
    ModelDataOutput_S2D_ME_s_RigidBody_Fric.h ModelDataOutput_S2D_ME_s_RigidBody_Fric.cpp
    TimeHistoryOutput_S2D_ME_s_RigidBody_Fric.h TimeHistoryOutput_S2D_ME_s_RigidBody_Fric.cpp
    # 2D CHM MPM
    Model_S2D_CHM_s.h Model_S2D_CHM_s.cpp
    Step_S2D_CHM_s.h
    Step_S2D_CHM_s.cpp
    Step_S2D_CHM_s_avg_stress.cpp
    Step_S2D_CHM_s_GIMP.cpp
    ModelDataOutput_S2D_CHM_s.h ModelDataOutput_S2D_CHM_s.cpp
    TimeHistoryOutput_S2D_CHM_s.h TimeHistoryOutput_S2D_CHM_s.cpp
    # 2D CHM FEM, u - U - p formualtion
    MatrixCoefficientSet.hpp
    Model_S2D_CHM_s_FEM_uUp.h Model_S2D_CHM_s_FEM_uUp.cpp
    Step_S2D_CHM_s_FEM_uUp.h Step_S2D_CHM_s_FEM_uUp.cpp
    TimeHistoryOutput_S2D_CHM_s_FEM_uUp.h
    )

add_subdirectory(TriangleMeshRigidBody/)

target_include_directories(
    SimulationCore
    PRIVATE
        ${SIMULATIONCORE_INCLUDE_DIR}
        ${UTILITIES_INCLUDE_DIR}
        ${EIGEN_INCLUDE_DIR}
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
