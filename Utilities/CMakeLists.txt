#===================================
#
# Utilities Project
#
#===================================
add_library(
    Utilities STATIC
    Utilities_pcp.h Utilities_pcp.cpp
    common_utilities.h common_utilities.cpp
    )

# Utilities headers
target_include_directories(
    Utilities PRIVATE
    ${UTILITIES_INCLUDE_DIR}
    # OpenGL
    ${GLAD_INCLUDE_DIR}
    ${GLFW_INCLUDE_DIR}
    ${GLM_INCLUDE_DIR}
    )

link_directories(${GLFW_LIBRARIES_DIR})
target_link_libraries(
    Utilities PRIVATE
    ${GLFW_LIBRARIES}
    )

add_subdirectory(GraphicUtils/)
add_subdirectory(MemoryUtils/)

if (MSVC)
    target_compile_definitions(Utilities PRIVATE _CRT_SECURE_NO_WARNINGS)
    target_compile_options(Utilities PRIVATE -YcUtilities_pcp.h)
    # exclude glad.c from precompilation
    # must add GraphicUtils/glad.c
    set_source_files_properties(GraphicUtils/glad.c PROPERTIES COMPILE_FLAGS /Y-)
endif()
