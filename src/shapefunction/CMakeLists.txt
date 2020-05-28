add_library(shapefunction STATIC "")

target_sources(shapefunction
    PRIVATE
        ${CMAKE_CURRENT_LIST_DIR}/shapefunction_fun.f90
        ${CMAKE_CURRENT_LIST_DIR}/shapefunction_obj.f90
        ${CMAKE_CURRENT_LIST_DIR}/shapefunction_wrapper.f90
    )

set_target_properties(shapefunction
    PROPERTIES
        POSITION_INDEPENDENT_CODE 1
    )

install(
    TARGETS
        shapefunction
    LIBRARY
        DESTINATION ${INSTALL_LIBDIR}
    OBJECTS
        DESTINATION ${INSTALL_OBJDIR}
    )
