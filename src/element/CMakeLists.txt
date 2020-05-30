add_library(element STATIC "")

target_sources(element
    PRIVATE
        ${CMAKE_CURRENT_LIST_DIR}/element2d_obj.f90
        ${CMAKE_CURRENT_LIST_DIR}/element3d_obj.f90
    )

set_target_properties(element
    PROPERTIES
        POSITION_INDEPENDENT_CODE 1
    )

target_link_libraries(element shapefunction integrand material)

install(
    TARGETS
        element
    LIBRARY
        DESTINATION ${INSTALL_LIBDIR}
    OBJECTS
        DESTINATION ${INSTALL_OBJDIR}
    )
