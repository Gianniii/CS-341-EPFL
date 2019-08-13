//=============================================================================
//
//   Exercise code for the lecture
//   "Introduction to Computer Graphics"
//   by Prof. Dr. Mario Botsch, Bielefeld University
//
//   Copyright (C) Computer Graphics Group, Bielefeld University.
//
//=============================================================================

//== INCLUDES =================================================================

#include "Cylinder.h"
#include "SolveQuadratic.h"

#include <array>
#include <cmath>

//== IMPLEMENTATION =========================================================

bool
Cylinder::
intersect(const Ray&  _ray,
          vec3&       _intersection_point,
          vec3&       _intersection_normal,
          double&     _intersection_t) const
{

   const vec3 projVtoAxis= dot(_ray.direction, axis) * axis;
   const vec3 diff = _ray.origin - center;
   const vec3 projDifftoAxis= dot(diff, axis)*axis;

   const double A = dot(_ray.direction - projVtoAxis, _ray.direction - projVtoAxis);
   const double B = 2 * dot(_ray.direction - projVtoAxis, diff - projDifftoAxis);
   const double C = dot((diff - projDifftoAxis), diff - projDifftoAxis) - radius*radius;

   std::array<double, 2> t;
   size_t sol = solveQuadratic(A,B,C,t);

   _intersection_t = NO_INTERSECTION;

    // Find the closest valid solution (in front of the viewer) && check if the point is in the cylinder (the projection of
    // the point in the axis is at most height/2 far away from the center)
    for (size_t i = 0; i < sol; ++i) {
        const vec3 projXtToAxis = dot(_ray(t[i]) - center,axis) * axis;
        const double dt = norm(projXtToAxis);
        if (t[i] > 0 && dt < height/2.0) _intersection_t = std::min(_intersection_t, t[i]);
    }

    if (_intersection_t == NO_INTERSECTION) return false;

    _intersection_point  = _ray(_intersection_t);

    const vec3 cPoint = _ray(_intersection_t) - center;
    const vec3 projXtToAxis = dot(cPoint,axis) * axis;
    double sgn = 1.0;
    if (sol >1) {
        if (_intersection_t == std::max(t[0], t[1])){
            sgn = -1.0;
        } 
    }
    
    _intersection_normal = sgn*normalize(cPoint - projXtToAxis);

    return true;
}
