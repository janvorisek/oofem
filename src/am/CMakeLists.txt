set (additiveman
    polygon_clipping.cpp
    gcode_program.cpp
    parser.cpp
    printer.cpp
    additivemanufacturingproblem.C
    freeconstantsurfaceload.C
    )

add_library (am
    ${LIB_TYPE}
    ${additiveman}
    )