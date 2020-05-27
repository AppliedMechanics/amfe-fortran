add_library(integrand STATIC "")

target_sources(integrand
    PRIVATE
        ${CMAKE_CURRENT_LIST_DIR}/integrand_structural2d_nonlinear.f90
        ${CMAKE_CURRENT_LIST_DIR}/integrand_structural3d_nonlinear.f90
        ${CMAKE_CURRENT_LIST_DIR}/fields.f90
    )

set_target_properties(integrand
    PROPERTIES
        POSITION_INDEPENDENT_CODE 1
    )

target_link_libraries(integrand tools)

install(
    TARGETS
        integrand
    LIBRARY
        DESTINATION ${INSTALL_LIBDIR}
    OBJECTS
        DESTINATION ${INSTALL_OBJDIR}
    )
