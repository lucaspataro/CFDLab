//
// Created by Moritz Gnisia on 05.04.20.
//

#ifndef CFDLAB_ENUMS_H
#define CFDLAB_ENUMS_H

enum class solver_type {
    NONE,
    SOR,
    MULTIGRID,
    CG,
    spLU
};

enum class velocity_type {
    U,
    V
};

enum class border_position {
    TOP,
    BOTTOM,
    LEFT,
    RIGHT
};

enum class neighbour_position {
    TOP,
    BOTTOM,
    LEFT,
    RIGHT
};

enum class cell_type {
    FLUID=1,
    NOSLIP,
    INLET,
    OUTLET,
    FREESLIP,
    LID
};

enum class boundary_type {
    B_N, 
    B_E, 
    B_S, 
    B_W, 
    B_NE, 
    B_NW, 
    B_SE, 
    B_SW,
    UNSET
};

enum class matrix_selection {
    ROW,
    COLUMN
};

enum class read_type {
    INT,
    DOUBLE
};

#endif //CFDLAB_ENUMS_H
