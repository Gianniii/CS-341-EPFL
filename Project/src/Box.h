#ifndef BOX_H
#define BOX_H

#include "vec3.h"
#include "Ray.h"

class Box {
public:
    /// Minimum point of the bounding box
    vec3 bb_min_;
    /// Maximum point of the bounding box
    vec3 bb_max_;
    /// longest axis of the bounding box
    int axis;
    Box(){
        bb_min_ = vec3(0, 0, 0);
        bb_max_ = vec3(0, 0, 0);
        axis=0;
    }

    Box(vec3 mn, vec3 mx) {
        bb_min_ = mn;
        bb_max_ = mx;
        getAxis();
    }

    /// merges two boxes and returns a new one
    Box merge(Box b1, Box b2) {
        return Box(min(b1.bb_min_, b2.bb_min_), max(b1.bb_max_, b2.bb_max_));
    }

    /// expands box bounds to fit the given new box as well
    void expand(Box bb) {
        bb_min_ = min(bb_min_, bb.bb_min_);
        bb_max_ = max(bb_max_, bb.bb_max_);
        getAxis();
    }

    /// Intersects box with ray
    bool intersect(const Ray &_ray) const
    {
        /**
        * Intersect the ray `_ray` with the axis-aligned bounding box of the mesh.
        * Note that the minimum and maximum point of the bounding box are stored
        * in the member variables `bb_min_` and `bb_max_`. Return whether the ray
        * intersects the bounding box.
        * This function is ued in `Mesh::intersect()` to avoid the intersection test
        * with all triangles of every mesh in the scene. The bounding boxes are computed
        * in `Mesh::compute_bounding_box()`.
        */
        //ray = o + t*d
        const vec3 &o = _ray.origin;
        const vec3 &d = _ray.direction;

        //min.x = o.x + t*d.x = bb_min.x * x_axis
        //min.y = o.y + t*d.y = bb_min.y * y_axis
        //min.z = o.z + t*d.z = bb_min.z * z_axis
        //similar for max coordinates

        double t_min_x = (bb_min_[0] - o[0]) / d[0];
        double t_max_x = (bb_max_[0] - o[0]) / d[0];

        //inspired by angleWeights function
        double t_min_min_x = std::min(t_min_x, t_max_x);
        double t_max_max_x = std::max(t_max_x, t_min_x);

        double t_min_y = (bb_min_[1] - o[1]) / d[1];
        double t_max_y = (bb_max_[1] - o[1]) / d[1];

        double t_min_min_y = std::min(t_min_y, t_max_y);
        double t_max_max_y = std::max(t_max_y, t_min_y);

        //ray misses box
        if ((t_min_min_x > t_max_max_y) || (t_min_min_y > t_max_max_x))
            return false;

        //take the one further away for the minimum intersection
        double t_min = std::max(t_min_min_x, t_min_min_y);
        //take the one closer for the maximum intersection
        double t_max = std::min(t_max_max_x, t_max_max_y);

        double t_min_z = (bb_min_[2] - o[2]) / d[2];
        double t_max_z = (bb_max_[2] - o[2]) / d[2];

        double t_min_min_z = std::min(t_min_z, t_max_z);
        double t_max_max_z = std::max(t_min_z, t_max_z);

        //misses box
        if ((t_min > t_max_max_z) || (t_min_min_z > t_max))
            return false;

        return true;
    }


private:
    void getAxis(){
        float xL = bb_max_[0] - bb_min_[0];
        float yL = bb_max_[1] - bb_min_[1];
        float zL = bb_max_[2] - bb_min_[2];
        float maxDist = fmax(xL, fmax(yL, zL));
        axis = maxDist == xL ? 0 : (maxDist == yL ? 1 : 2);
    }

};

#endif