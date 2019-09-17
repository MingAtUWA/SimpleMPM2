#===================================
#
# Utilities/GraphicUtils
#
#===================================
target_sources(
    Utilities PRIVATE
    # include glad, glfw, glm
    OpenGL_headers.h
    # glad
    glad.c
    # opengl shader program
    ShaderProgram.h ShaderProgram.cpp
    # gif
    gif.h
    )
