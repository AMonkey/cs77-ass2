#include "primitive.h"

///@file igl/primitive.cpp Primitives. @ingroup igl

/// Computes transformation matrix of TransformedSurface
/// @param transformed The given TransformedSurface for which to compute a transformation matrix
/// @return The transformation matrix
mat4f transformed_matrix(TransformedSurface* transformed) {
    auto m = identity_mat4f;

    put_your_code_here("Transformation");

    return m;
}

/// Computes inverse transformation matrix of TransformedSurface
/// @param transformed The given TransformedSurface for which to compute an inverse transformation matrix
/// @return The inverse transformation matrix
mat4f transformed_matrix_inv(TransformedSurface* transformed) {
    auto mi = identity_mat4f;

    put_your_code_here("Transformation (inverse)");
    
    return mi;
}

