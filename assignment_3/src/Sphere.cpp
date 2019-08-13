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

#include "Sphere.h"
#include "SolveQuadratic.h"

//== IMPLEMENTATION =========================================================


Sphere::Sphere(const vec3& _center, double _radius)
: center(_center), radius(_radius)
{
}


//-----------------------------------------------------------------------------


bool
Sphere::
intersect(const Ray& ray, vec3& intersection_point, vec3& intersection_normal, double& intersection_t) const
{
    /*
    We are solving the equation:
        t^2 * ||d||^2 + t * (2d dot (o - c)) + ||o-c||^2 - r^2 = 0
    for the intersection distance "t", where
        c = center
        d = ray.direction
        o = ray.origin
        r = radius
    */
    const vec3 &d = ray.direction;
    const vec3 oc = ray.origin - center;
    std::array<double, 2> t;

    size_t number_of_solutions = solveQuadratic(
        dot(d, d),                     // t^2 * ||d||^2
        2 * dot(d, oc),                // t * (2d dot (o - c))
        dot(oc, oc) - radius * radius, // ||o-c||^2 - r^2
        t                              // where to store solutions
    );

    // Initialize the intersection distance to "infinity" (indicating no intersection)
    intersection_t = NO_INTERSECTION;

    // Find the closest valid (in front of the viewer) intersection, if one exists
    for (size_t i = 0; i < number_of_solutions; i++) {
        const double t_solution = t[i];

        // t <= 0 means the intersection is behind the camera
        if ((t_solution > 0) && (t_solution < intersection_t)) {
            intersection_t = t_solution;
        }
    }

    // Nothing to do if we don't have an intersection...
    if (intersection_t == NO_INTERSECTION) return false;

    // Otherwise, calculate information about the intersection
    intersection_point  = ray(intersection_t);
    intersection_normal = (intersection_point - center) / radius;

    return true;
}

//=============================================================================
