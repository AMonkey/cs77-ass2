#include "primitive.h"

///@file igl/primitive.cpp Primitives. @ingroup igl

/// Computes transformation matrix of TransformedSurface
/// @param transformed The given TransformedSurface for which to compute a transformation matrix
/// @return The transformation matrix
mat4f transformed_matrix(TransformedSurface* transformed) {
    auto t = transformed->translation;
    auto r = transformed->rotation_euler;
    auto s = transformed->scale;

    // Translation
    auto mat_t = mat4f(vec4f(1, 0, 0, t.x),
                       vec4f(0, 1, 0, t.y),
                       vec4f(0, 0, 1, t.z),
                       vec4f(0, 0, 0, 1)

                 );

    // Scale
    auto mat_s = mat4f(vec4f(s.x, 0, 0, 0),
                       vec4f(0, s.y, 0, 0),
                       vec4f(0, 0, s.z, 0),
                       vec4f(0, 0, 0, 1)

                 );

    // Rotation about z
    auto mat_r_z = mat4f(vec4f(cos(r.z), -1*sin(r.z), 0, 0),
                         vec4f(sin(r.z), cos(r.z), 0, 0),
                         vec4f(0, 0, 1, 0),
                         vec4f(0, 0, 0, 1)

                   );

    // Rotation about y
    auto mat_r_y = mat4f(vec4f(cos(r.y), 0, sin(r.y), 0),
                         vec4f(0, 1, 0, 0),
                         vec4f(-1*sin(r.y), 0, cos(r.y), 0),
                         vec4f(0, 0, 0, 1)

                   );

    // Rotation about x
    auto mat_r_x = mat4f(vec4f(1, 0, 0, 0),
                         vec4f(0, cos(r.x), -1*sin(r.x), 0),
                         vec4f(0, sin(r.x), cos(r.x), 0),
                         vec4f(0, 0, 0, 1)

                   );

    auto m = mat_t * mat_r_z * mat_r_y * mat_r_x * mat_s;

    return m;
    //return m * transformed_matrix_inv(transformed); // uncomment to check inv

}

/// Computes inverse transformation matrix of TransformedSurface
/// @param transformed The given TransformedSurface for which to compute an inverse transformation matrix
/// @return The inverse transformation matrix
mat4f transformed_matrix_inv(TransformedSurface* transformed) {
    auto t = transformed->translation;
    auto r = transformed->rotation_euler;
    auto s = transformed->scale;

    // To invert, subtract instead of adding
    auto mat_t = mat4f(vec4f(1, 0, 0, -1*t.x),
                       vec4f(0, 1, 0, -1*t.y),
                       vec4f(0, 0, 1, -1*t.z),
                       vec4f(0, 0, 0, 1)

                 );

    // To invert, take ^-1 of scale factors
    auto mat_s = mat4f(vec4f((1/s.x), 0, 0, 0),
                       vec4f(0, (1/s.y), 0, 0),
                       vec4f(0, 0, (1/s.z), 0),
                       vec4f(0, 0, 0, 1)

                 );

    // To invert, take opposite of theta (only affects sin as cos is even)
    auto mat_r_z = mat4f(vec4f(cos(r.z), sin(r.z), 0, 0),
                         vec4f(-1*sin(r.z), cos(r.z), 0, 0),
                         vec4f(0, 0, 1, 0),
                         vec4f(0, 0, 0, 1)

                   );

    auto mat_r_y = mat4f(vec4f(cos(r.y), 0, -1*sin(r.y), 0),
                         vec4f(0, 1, 0, 0),
                         vec4f(sin(r.y), 0, cos(r.y), 0),
                         vec4f(0, 0, 0, 1)

                   );

    auto mat_r_x = mat4f(vec4f(1, 0, 0, 0),
                         vec4f(0, cos(r.x), sin(r.x), 0),
                         vec4f(0, -1*sin(r.x), cos(r.x), 0),
                         vec4f(0, 0, 0, 1)

                   );
    return mat_s * mat_r_x * mat_r_y * mat_r_z * mat_t;

}

