#ifndef TRIANGLE_H
#define TRIANGLE_H

#include "vec3.h"
#include "Box.h"

/// a triangle is specified by three indices and a normal
struct Triangle
{
    /// index of first vertex (for array Mesh::vertices_)
    int i0;
    /// index of second vertex (for array Mesh::vertices_)
    int i1;
    /// index of third vertex (for array Mesh::vertices_)
    int i2;
    /// triangle normal
    vec3 normal;
    /// triangle midpoint
    vec3 midpoint;
    /// triangle's bounding box
    Box bounds;
};

#endif