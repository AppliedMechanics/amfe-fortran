add_library(quadrature STATIC "")

target_sources(quadrature
    PRIVATE
        ${CMAKE_CURRENT_LIST_DIR}/quadrature.f90
    )

set_target_properties(quadrature
    PROPERTIES
        POSITION_INDEPENDENT_CODE 1
    )

install(
    TARGETS
        quadrature
    LIBRARY
        DESTINATION ${INSTALL_LIBDIR}
    OBJECTS
        DESTINATION ${INSTALL_OBJDIR}
    )
