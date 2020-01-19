#===================================
#
# ConstitutiveModel Project
#
#===================================
add_library(
    ConstitutiveModel STATIC
    ConstitutiveModel_pcp.h ConstitutiveModel_pcp.cpp
    ConstitutiveModel.h
    # Linear Elasticity
    LinearElasticity.h LinearElasticity.cpp
    # Modified Cam Clay
    ModifiedCamClay.h ModifiedCamClay.cpp
    )

if (MSVC)
    target_compile_definitions(ConstitutiveModel PRIVATE _CRT_SECURE_NO_WARNINGS)
    target_compile_options(ConstitutiveModel PRIVATE -YcConstitutiveModel_pcp.h)
endif()
