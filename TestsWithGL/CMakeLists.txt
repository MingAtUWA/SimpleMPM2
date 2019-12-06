#===================================
#
# TestsWithGL Project
#
#===================================
add_executable(
    TestsWithGL TestsWithGL_main.cpp
    TestsWithGL_pcp.h TestsWithGL_pcp.cpp
    test_sim_core.h
    test_post_processor.h
    # Test functions
    # test Step solve
    test_solve_functions.cpp
    # rigid body testing
    DrawTriangleMesh.h DrawTriangleMesh.cpp
    display_triangle_mesh.cpp
    test_triangle_mesh_circle.cpp
    test_triangle_mesh_square.cpp
    test_rigid_body_square.cpp
    # mpm - rigid body contact
    test_mpm_rigidbody_circle.cpp
    test_mpm_rigidbody_square.cpp
    test_mpm_rigidbody_cantilever.cpp
    test_mpm_rigidbody_cantilever_fric.cpp
    test_mpm_rigidbody_bar_compression.cpp
    test_slide_down_frictional_slope.cpp
    test_slide_down_frictional_slope2.cpp
    # chm - mpm
    test_mpm_chm_s_1d_consolidation.cpp
    test_t2d_mpm_chm_s_1d_consolidation.cpp
    test_t2d_mpm_chm_s_t_bar.cpp
    # Animation output
    test_mpm_rigidbody_animation.cpp
    # Generate material points from mesh
    test_get_pcls_from_mesh.cpp
    # GIMP
    test_init_pcl_gimp.cpp
    # implicit FEM or MPM
    test_matrix_coefficient_set.cpp
    test_fem_me_s_up_1dbar.cpp
    test_fem_chm_s_1d_consolidation.cpp
    test_imp_mpm_me_s_up_1dbar.cpp
    test_imp_mpm_chm_s_uup_1d_consolidation.cpp
    # ColorGraph
    test_color_graph.cpp
    # Search triangle
    test_triangle_searching.cpp
    )

target_include_directories(
    TestsWithGL
    PRIVATE
    ${SIMULATIONCORE_INCLUDE_DIR}
    ${CMAKE_SOURCE_DIR}/PostProcessor
    ${UTILITIES_INCLUDE_DIR}
    ${EIGEN_INCLUDE_DIR}
    ${OPENGL_INCLUDE_DIR}
    )

link_directories(${GLFW_LIBRARIES_DIR})
target_link_libraries(
    TestsWithGL
    PRIVATE
        SimulationCore
        PostProcessor
        Utilities
        ${GLFW_LIBRARIES}
    )

if (WIN32)
    # supress warning for certain functions in visual studio
    target_compile_definitions(TestsWithGL PRIVATE _CRT_SECURE_NO_WARNINGS)
    # use precompiled header
    target_compile_options(TestsWithGL PRIVATE -YcTestsWithGL_pcp.h)
endif()
