add_library(shapefunction SHARED "")

target_sources(shapefunction
    PRIVATE
        ${CMAKE_CURRENT_LIST_DIR}/shapefunction_fun.f90
        ${CMAKE_CURRENT_LIST_DIR}/shapefunction_wrapper.f90
    )

set_target_properties(shapefunction
    PROPERTIES
        POSITION_INDEPENDENT_CODE 1
    )

install(
    TARGETS
        shapefunction
    ARCHIVE
        DESTINATION ${INSTALL_LIBDIR}
    RUNTIME
        DESTINATION ${INSTALL_LIBDIR}
    LIBRARY
        DESTINATION ${INSTALL_LIBDIR}
    )
