add_library(material STATIC "")

target_sources(material
    PRIVATE
        ${CMAKE_CURRENT_LIST_DIR}/material_fun.f90
        ${CMAKE_CURRENT_LIST_DIR}/material_obj.f90
    )

set_target_properties(material
    PROPERTIES
        POSITION_INDEPENDENT_CODE 1
    )

target_link_libraries(material tools)

install(
    TARGETS
        material
    LIBRARY
        DESTINATION ${INSTALL_LIBDIR}
    OBJECTS
        DESTINATION ${INSTALL_OBJDIR}
    )
