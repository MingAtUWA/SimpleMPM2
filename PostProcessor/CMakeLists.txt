# ===================================
#
# PostProcessor project
#
# ===================================
add_library(
    PostProcessor STATIC 
    PostProcessor_pcp.h PostProcessor_pcp.cpp
    GenerateAnimation.h GenerateAnimation.cpp
    GA_S2D_ME_s_RigidBody.h GA_S2D_ME_s_RigidBody.cpp
    GA_S2D_ME_s_RigidBody_Fric.h GA_S2D_ME_s_RigidBody_Fric.cpp
    GA_S2D_CHM_s.h GA_S2D_CHM_s.cpp
    GA_S2D_ME_s_up.h GA_S2D_ME_s_up.cpp
    GA_S2D_CHM_s_uUp.h GA_S2D_CHM_s_uUp.cpp
    GA_T2D_CHM_s.h GA_T2D_CHM_s.cpp
    DisplayModel.h DisplayModel.cpp
    DisplayModel_T2D.h DisplayModel_T2D.cpp
    SquareParticleSystem.h SquareParticleSystem.cpp
    # color particle according to variables
    ColorGraph.h
    ColorfulSquareParticleSystem.h ColorfulSquareParticleSystem.cpp
    GA_T2D_CHM_s_color.h GA_T2D_CHM_s_color.cpp
    )

target_include_directories(
    PostProcessor PRIVATE
    ${UTILITIES_INCLUDE_DIR}
    ${SIMULATIONCORE_INCLUDE_DIR}
    ${OPENGL_INCLUDE_DIR}
    )

target_link_libraries(
    PostProcessor PRIVATE
    Utilities
    )

if (MSVC)
    target_compile_definitions(PostProcessor PRIVATE _CRT_SECURE_NO_WARNINGS)
    target_compile_options(PostProcessor PRIVATE -YcPostProcessor_pcp.h)
endif()
