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
    ResultFile.h ResultFileType.h
    # Plain binary result file
    ResultFile_PlainBin.h ResultFile_PlainBin_DataStruct.h
    # XML result file
    ResultFile_XML.h
    # hdf5 result file
    ResultFile_hdf5.h ResultFile_hdf5.cpp
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
    # Implicit time integration
    MatrixCoefficientSet.hpp
    # 2D ME FEM u - p formulation
    Model_S2D_ME_s_FEM_up.h Model_S2D_ME_s_FEM_up.cpp
    Step_S2D_ME_s_FEM_up.h Step_S2D_ME_s_FEM_up.cpp
    ModelDataOutput_S2D_ME_s_FEM_up.h ModelDataOutput_S2D_ME_s_FEM_up.cpp
    TimeHistoryOutput_S2D_ME_s_FEM_up.h TimeHistoryOutput_S2D_ME_s_FEM_up.cpp
    # 2D CHM FEM, u - U - p formulation
    Model_S2D_CHM_s_FEM_uUp.h Model_S2D_CHM_s_FEM_uUp.cpp
    Step_S2D_CHM_s_FEM_uUp.h Step_S2D_CHM_s_FEM_uUp.cpp Step_S2D_CHM_s_FEM_uUp_mat_and_vec.cpp
    ModelDataOutput_S2D_CHM_s_FEM_uUp.h ModelDataOutput_S2D_CHM_s_FEM_uUp.cpp
    TimeHistoryOutput_S2D_CHM_s_FEM_uUp.h TimeHistoryOutput_S2D_CHM_s_FEM_uUp.cpp
    # 2D ME MPM, u - p formulation
    Model_S2D_ME_s_up.h Model_S2D_ME_s_up.cpp
    Step_S2D_ME_s_up.h Step_S2D_ME_s_up.cpp
    Step_S2D_ME_s_up_reequilibration.cpp Step_S2D_ME_s_up_mat_and_vec.cpp
    ModelDataOutput_S2D_ME_s_up.h ModelDataOutput_S2D_ME_s_up.cpp
    TimeHistoryOutput_S2D_ME_s_up.h TimeHistoryOutput_S2D_ME_s_up.cpp
    # 2D CHM MPM, u - U - p formulation
    Model_S2D_CHM_s_uUp.h Model_S2D_CHM_s_uUp.cpp
    Step_S2D_CHM_s_uUp.h Step_S2D_CHM_s_uUp.cpp
    Step_S2D_CHM_s_uUp_form_and_solve.cpp Step_S2D_CHM_s_uUp_mat_and_vec.cpp
    ModelDataOutput_S2D_CHM_s_uUp.h ModelDataOutput_S2D_CHM_s_uUp.cpp
    TimeHistoryOutput_S2D_CHM_s_uUp.h TimeHistoryOutput_S2D_CHM_s_uUp.cpp
    # Triangular 2D Mesh
    # 2D ME MPM
    Model_T2D_ME_s.h Model_T2D_ME_s.cpp
    Step_T2D_ME_s.h Step_T2D_ME_s.cpp
    Step_T2D_ME_s_Geostatic.h Step_T2D_ME_s_Geostatic.cpp
    #ModelDataOutput_T2D_ME_s.h ModelDataOutput_T2D_ME_s.cpp
    #TimeHistoryOutput_T2D_ME_s.h TimeHistoryOutput_T2D_ME_s.cpp
    # 2D CHM MPM
    Model_T2D_CHM_s.h Model_T2D_CHM_s.cpp
    Step_T2D_CHM_s.h Step_T2D_CHM_s.cpp
    ModelDataOutput_T2D_CHM_s.h ModelDataOutput_T2D_CHM_s.cpp
    TimeHistoryOutput_T2D_CHM_s.h TimeHistoryOutput_T2D_CHM_s.cpp
    # Rigid ciricle controlled by displacement
    DispConRigidCircle.h DispConRigidCircle.cpp
    # Strain enhancement
    Step_T2D_CHM_s_SE.h Step_T2D_CHM_s_SE.cpp
    Step_T2D_CHM_s_SE_Geostatic.h Step_T2D_CHM_s_SE_Geostatic.cpp
    TimeHistoryOutput_T2D_CHM_s_SE.h TimeHistoryOutput_T2D_CHM_s_SE.cpp
    ModelDataOutput_T2D_ME_s.h ModelDataOutput_T2D_ME_s.cpp
    TimeHistoryOutput_T2D_ME_s.h TimeHistoryOutput_T2D_ME_s.cpp
    TimeHistoryOutput_T2D_CHM_s_SE_Geostatic.h TimeHistoryOutput_T2D_CHM_s_SE_Geostatic.cpp
    # output mc container to hdf5 file
    output_model_container_to_hdf5.h output_model_container_to_hdf5.cpp
    )

add_subdirectory(TriangleMeshRigidBody/)

target_include_directories(
    SimulationCore PRIVATE
    ${SIMULATIONCORE_INCLUDE_DIR}
    ${CONSTITUTIVEMODEL_INCLUDE_DIR}
    ${UTILITIES_INCLUDE_DIR}
    ${EIGEN_INCLUDE_DIR}
    ${OPENGL_INCLUDE_DIR}
    ${HDF5_INCLUDE_DIR}
    )

target_link_libraries(
    SimulationCore PRIVATE
    Utilities
    ConstitutiveModel
    ${HDF5_LIBRARIES}
    )

if (MSVC)
    target_compile_definitions(SimulationCore PRIVATE _CRT_SECURE_NO_WARNINGS)
    target_compile_options(SimulationCore PRIVATE -YcSimulationCore_pcp.h)
endif()
