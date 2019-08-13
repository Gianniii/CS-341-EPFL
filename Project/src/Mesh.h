//=============================================================================
//
//   Exercise code for the lecture
//   "Introduction to Computer Graphics"
//   by Prof. Dr. Mario Botsch, Bielefeld University
//
//   Copyright (C) Computer Graphics Group, Bielefeld University.
//
//=============================================================================

#ifndef MESH_H
#define MESH_H


//== INCLUDES =================================================================

#include "Object.h"
#include <vector>
#include <utility>
#include <string>
#include "Triangle.h"
#include "Box.h"

//== CLASS DEFINITION =========================================================


/// \class Mesh Mesh.h
/// This class represents a simple triangle mesh, stored as an indexed face set,
/// i.e., as an array of vertices and an array of triangles.
class KDNode;
class Mesh : public Object
{
public:
    KDNode* root;

    /// Array of triangles
    std::vector<Triangle> triangles_;

    /// This type is used to choose between flat shading and Phong shading
    enum Draw_mode {FLAT, PHONG};

    /// Construct a mesh by parsing its path and properties from an input
    /// stream. The mesh path read from the file is relative to the 
    /// scene file's path "scenePath".
    Mesh(std::istream &is, const std::string &scenePath);

    /// Intersect mesh with ray (calls ray-triangle intersection)
    /// If \c _ray intersects a face of the mesh, it provides the following results:
    /// \param[in] _ray the ray to intersect the mesh with
    /// \param[out] _intersection_point the point of intersection
    /// \param[out] _intersection_normal the surface normal at intersection point
    /// \param[out] _intersection_t ray parameter at the intersection point
    virtual bool intersect(const Ray& _ray,
                           vec3&      _intersection_point,
                           vec3&      _intersection_normal,
                           double&    _intersection_t,
                           vec3& tex) const override;

private:
    /// a vertex consists of a position and a normal
    struct Vertex
    {
        /// vertex position
        vec3 position;
        /// vertex normal
        vec3 normal;
        vec3 textures;
    };

public:
    /// Read mesh from an OFF file
    bool read(const std::string &_filename);

    /// Compute normal vectors for triangles and vertices
    void compute_normals();

    /// Compute the axis-aligned bounding box, store minimum and maximum point in bb_min_ and bb_max_
    void compute_bounding_box();

    /// Does \c _ray intersect the bounding box of the mesh?
    bool intersect_bounding_box(const Ray& _ray) const;

    /// Intersect a triangle with a ray. Return whether there is an intersection.
    /// If there is an intersection, store intersection data.
    /// This function overrides Object::intersect().
    /// \param[in] _triangle the triangle to be intersected
    /// \param[in] _ray the ray to intersect the triangle with
    /// \param[out] _intersection_point the point of intersection
    /// \param[out] _intersection_normal the surface normal at intersection point
    /// \param[out] _intersection_t ray parameter at the intersection point
    bool intersect_triangle(const Triangle&  _triangle,
                            const Ray&       _ray,
                            vec3&            _intersection_point,
                            vec3&            _intersection_normal,
                            double&          _intersection_t,
                            vec3& tex) const;

private:
    /// Does this mesh use flat or Phong shading?
    Draw_mode draw_mode_;

    /// Array of vertices
    std::vector<Vertex> vertices_;

    /// Minimum point of the bounding box
    vec3 bb_min_;
    /// Maximum point of the bounding box
    vec3 bb_max_;
};

class KDNode
{
public:
    /// node's bounding box
    Box bbox;
    /// left child node
    KDNode* left;
    /// right child node
    KDNode* right;
    /// indices of triangles in this box (only filled if leaf node)
    std::vector<int> triangles;
    /// indicator whether it's a leaf node
    bool leaf;

    /// returns if node is empty
    bool isEmpty() const{
        return triangles.size() == 0;
    }

    // builds a new node referencing to the triangles of the passed mesh
    static KDNode* build(Mesh* m, std::vector<int>& tris, int depth)
    {
        KDNode* node = new KDNode();
        
        if (tris.size() == 0) {
            node->triangles = tris;
            node->leaf = true;
            return node;
        }
        node->bbox = m->triangles_[tris[0]].bounds;

        if(tris.size() == 1) {
            node->leaf = true;
            node->triangles = tris;
            node->left = NULL;
            node->right = NULL;
            return node;
        }

        vec3 midpoint = vec3(0, 0, 0);
        for (int i: tris) {
            node->bbox.expand(m->triangles_[i].bounds);
            midpoint = midpoint + (m->triangles_[i].midpoint * 1.0/tris.size());
        }

        std::vector<int> left_tris;
        std::vector<int> right_tris;
        left_tris.reserve(tris.size());
        right_tris.reserve(tris.size());
        
        for (int i: tris){
            switch (node->bbox.axis) {
                case 0:
                    midpoint[0] >= m->triangles_[i].midpoint[0] ? right_tris.push_back(i) : left_tris.push_back(i);
                    break;
                case 1:
                    midpoint[1] >= m->triangles_[i].midpoint[1] ? right_tris.push_back(i) : left_tris.push_back(i);
                    break;
                case 2:
                    midpoint[2] >= m->triangles_[i].midpoint[2] ? right_tris.push_back(i) : left_tris.push_back(i);
                    break;
                default:
                    std::cerr << "axis is " << node->bbox.axis;
                    break;
            }
        }

        if(depth < 10/*(float)matches / left_tris.size() < 0.5 && (float)matches / right_tris.size() < 0.5*/){
            node->left = build(m, left_tris, depth+1);
            node->right = build(m, right_tris, depth+1);
        } else {
            node->triangles = tris;
            node->leaf = true;
            node->left = NULL;
            node->right = NULL;
        }

        return node;
    }

    /// check if ray intersects with bounding box of the node and return pointer to index vector if leaf node
    bool hit(const Ray &_ray, std::vector<std::pair<int*,int>> &tris)
    {
        if(bbox.intersect(_ray)){
		    bool hit_tri = false;
            if (leaf) {
                if (isEmpty()) {
                    return false;
                } else {
                    tris.push_back(std::make_pair(triangles.data(),triangles.size()));
                    return true;
                }
            }
            else
            {
                bool hitLeft = left->hit(_ray, tris);
                bool hitRight = right->hit(_ray, tris);
                return hitLeft || hitRight;
            }
        } else {
            return false;
        }
    }
};

//=============================================================================
#endif // MESH_H defined
//=============================================================================
