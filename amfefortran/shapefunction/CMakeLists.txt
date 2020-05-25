add_library(shapefunction "")

target_sources(shapefunction
    PUBLIC
        ${CMAKE_CURRENT_LIST_DIR}/shapefunction_fun.f90
    )
