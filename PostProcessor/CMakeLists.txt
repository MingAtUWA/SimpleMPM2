# ===================================
#
# PostProcessor project
#
# ===================================
add_library(
    PostProcessor STATIC 
    PostProcessor_pcp.h PostProcessor_pcp.cpp
    GenerateAnimation.h GenerateAnimation.cpp
    )

target_include_directories(
    PostProcessor PRIVATE
    ${UTILITIES_INCLUDE_DIR}
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
