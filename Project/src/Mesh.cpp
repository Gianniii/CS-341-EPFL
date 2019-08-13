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

#include "Mesh.h"
#include <fstream>
#include <string>
#include <stdexcept>
#include <limits>

#include "tiny_obj_loader.h"

//== IMPLEMENTATION ===========================================================

Mesh::Mesh(std::istream &is, const std::string &scenePath)
{
	std::string meshFile, mode;
	is >> meshFile;

	// load mesh from file
	read(scenePath.substr(0, scenePath.find_last_of("/\\") + 1) + meshFile); // Use both Unix and Windows path separators

	is >> mode;
	if (mode == "FLAT")
		draw_mode_ = FLAT;
	else if (mode == "PHONG")
		draw_mode_ = PHONG;
	else
		throw std::runtime_error("Invalid draw mode " + mode);

	is >> material;
}

//-----------------------------------------------------------------------------

bool Mesh::read(const std::string &_filename)
{
	// read a mesh in OBJ format
	if (_filename.substr(_filename.find_last_of(".") + 1) == "obj")
	{

		tinyobj::attrib_t attrib;
		std::vector<tinyobj::shape_t> shapes;
		std::vector<tinyobj::material_t> materials;

		std::string warn;
		std::string err;

		// load an Obj file from _filename path
		bool ret = tinyobj::LoadObj(&attrib, &shapes, &materials, &warn, &err, _filename.c_str());


		std::string mattt = materials[0].diffuse_texname;

		// test if the obj file has a diffuse texture, if so, create a Texture object and store it in the texture optional parameter
		// of the Material of the object
		if (!materials[0].diffuse_texname.empty())
		{

			material.texture.reset(new Texture(materials[0].diffuse_texname));
		}

		if (!warn.empty())
		{
			std::cout << warn << std::endl;
		}

		if (!err.empty())
		{
			std::cerr << err << std::endl;
		}

		if (!ret)
		{
			exit(1);
		}

		// laod vertices from the obj file and store them in the vertices array of the Mesh Object
		// load textures and position
		Vertex vert;
		vertices_.clear();

		for (size_t v = 0; v < attrib.vertices.size() / 3; v++)
		{
			tinyobj::real_t vx = attrib.vertices[3 * v + 0];
			tinyobj::real_t vy = attrib.vertices[3 * v + 1];
			tinyobj::real_t vz = attrib.vertices[3 * v + 2];
			tinyobj::real_t tx = attrib.texcoords[2*v+0];
      		tinyobj::real_t ty = attrib.texcoords[2*v+1];
			vert.position = vec3(vx, vy, vz);
			vert.textures = vec3(tx,ty,0);
			vertices_.push_back(vert);
		};

		// Store indices of vertices of each triangle of the Mesh
		Triangle t;
		triangles_.clear();
		triangles_.reserve(shapes[0].mesh.num_face_vertices.size());

		for (size_t s = 0; s < shapes.size(); s++)
		{
			 
			// Loop over faces(polygon)
			size_t index_offset = 0;
			size_t tri;
				for (size_t f = 0; f < shapes[s].mesh.num_face_vertices.size(); f++)
			{
				int fv = shapes[s].mesh.num_face_vertices[f];
				// Loop over vertices in the face.
				for (size_t v = 0; v < fv; v++)
				{
					tinyobj::index_t idx = shapes[s].mesh.indices[index_offset + v];
					switch (v)
					{
					case 0:
						t.i0 = idx.vertex_index;
						break;
					case 1:

						t.i1 = idx.vertex_index;
						break;
					case 2:

						t.i2 = idx.vertex_index;
						break;

					default:
						break;
					}

				}
				index_offset += fv;
				triangles_.push_back(t);
			}
		}

		
	}
	//else read off file
	else
	{
		std::ifstream ifs(_filename);
		if (!ifs)
		{
			std::cerr << "Can't open " << _filename << "\n";
			return false;
		}

		// read OFF header
		std::string s;
		unsigned int nV, nF, dummy, i;
		ifs >> s;
		if (s != "OFF")
		{
			std::cerr << "No OFF file\n";
			return false;
		}
		ifs >> nV >> nF >> dummy;
		std::cout << "\n  read " << _filename << ": " << nV << " vertices, " << nF << " triangles";

		// read vertices
		Vertex v;
		vertices_.clear();
		vertices_.reserve(nV);
		for (i = 0; i < nV; ++i)
		{
			ifs >> v.position;


			vertices_.push_back(v);
		}

		// read triangles
		Triangle t;
		triangles_.clear();
		triangles_.reserve(nF);
		for (i = 0; i < nF; ++i)
		{
			ifs >> dummy >> t.i0 >> t.i1 >> t.i2;
			triangles_.push_back(t);
		}

		// close file
		ifs.close();
	}
	std::vector<int> tris;
	tris.reserve(triangles_.size());
	//Building KDTree for Mesh
	for(int i = 0; i < triangles_.size(); ++i) {
	 	tris.push_back(i);
	}

	compute_normals();

	compute_bounding_box();

	root = KDNode::build(this, tris, 0);

	return true;
}

//-----------------------------------------------------------------------------

// Determine the weights by which to scale triangle (p0, p1, p2)'s normal when
// accumulating the vertex normals for vertices 0, 1, and 2.
// (Recall, vertex normals are a weighted average of their incident triangles'
// normals, and in our raytracer we'll use the incident angles as weights.)
// \param[in] p0, p1, p2    triangle vertex positions
// \param[out] w0, w1, w2    weights to be used for vertices 0, 1, and 2
void angleWeights(const vec3 &p0, const vec3 &p1, const vec3 &p2,
				  double &w0, double &w1, double &w2)
{
	// compute angle weights
	const vec3 e01 = normalize(p1 - p0);
	const vec3 e12 = normalize(p2 - p1);
	const vec3 e20 = normalize(p0 - p2);
	w0 = acos(std::max(-1.0, std::min(1.0, dot(e01, -e20))));
	w1 = acos(std::max(-1.0, std::min(1.0, dot(e12, -e01))));
	w2 = acos(std::max(-1.0, std::min(1.0, dot(e20, -e12))));
}

//-----------------------------------------------------------------------------

void Mesh::compute_normals()
{
	for (Vertex &v : vertices_)
	{
		v.normal = vec3(0, 0, 0);
	}

	// compute triangle normals
	for (Triangle &t : triangles_)
	{
		const vec3 &p0 = vertices_[t.i0].position;
		const vec3 &p1 = vertices_[t.i1].position;
		const vec3 &p2 = vertices_[t.i2].position;
		t.normal = normalize(cross(p1 - p0, p2 - p0));

		float xmin = fmin(p0[0], fmin(p1[0], p2[0]));
        float ymin = fmin(p0[1], fmin(p1[1], p2[1]));
        float zmin = fmin(p0[2], fmin(p1[2], p2[2]));

        float xmax = fmax(p0[0], fmax(p1[0], p2[0]));
        float ymax = fmax(p0[1], fmax(p1[1], p2[1]));
        float zmax = fmax(p0[2], fmax(p1[2], p2[2]));

		t.bounds = Box(vec3(xmin, ymin, zmin), vec3(xmax, ymax, zmax));

		t.midpoint = (p0 + p1 + p2) / 3;

		/**
         * In some scenes (e.g the office scene) some objects should be flat
         * shaded (e.g. the desk) while other objects should be Phong shaded to appear
         * realistic (e.g. chairs). You have to implement the following:
         * - Compute vertex normals by averaging the normals of their incident triangles.
         * - Store the vertex normals in the Vertex::normal member variable.
         * - Weigh the normals by their triangles' angles.
         */

		double w0 = 0, w1 = 0, w2 = 0;
		angleWeights(p0, p1, p2, w0, w1, w2);

		vertices_[t.i0].normal += w0 * t.normal;
		vertices_[t.i1].normal += w1 * t.normal;
		vertices_[t.i2].normal += w2 * t.normal;
	}

	// normalize vector normals
	for (Vertex &v : vertices_)
	{
		v.normal = normalize(v.normal);
	}
}

//-----------------------------------------------------------------------------

void Mesh::compute_bounding_box()
{
	bb_min_ = vec3(std::numeric_limits<double>::max());
	bb_max_ = vec3(std::numeric_limits<double>::lowest());

	for (Vertex v : vertices_)
	{
		bb_min_ = min(bb_min_, v.position);
		bb_max_ = max(bb_max_, v.position);
	}
}

//-----------------------------------------------------------------------------

bool Mesh::intersect_bounding_box(const Ray &_ray) const
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

//-----------------------------------------------------------------------------

bool Mesh::intersect(const Ray &_ray,
					 vec3 &_intersection_point,
					 vec3 &_intersection_normal,
					 double &_intersection_t,
					 vec3 &tex) const
{
	// check bounding box intersection
	// if (!intersect_bounding_box(_ray))
	// {
	// 	return false;
	// }
	std::vector<std::pair<int*,int>> tris;
	tris.reserve(triangles_.size());
	if (!root->hit(_ray, tris)) {
		return false;
	}

	vec3 p, n;
	double t;

	_intersection_t = NO_INTERSECTION;

	// for each triangle
	// for (const Triangle &triangle : triangles_)
	// {
	for (std::pair<int*, int> pair : tris) {
		for (unsigned i = 0; i < pair.second; ++i){
			Triangle triangle = triangles_[*pair.first++];
			// does ray intersect triangle?
			if (intersect_triangle(triangle, _ray, p, n, t, tex))
			{
				// is intersection closer than previous intersections?
				if (t < _intersection_t)
				{
					// store data of this intersection
					_intersection_t = t;
					_intersection_point = p;
					_intersection_normal = n;
				}
			}
		}
	}

	return (_intersection_t != NO_INTERSECTION);
}

//-----------------------------------------------------------------------------

double determ3x3(const vec3 &a, const vec3 &b, const vec3 &c)
{
	//return a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0]) + a[2] * (b[0] * c[1] - b[1] * c[0]);
	vec3 p = cross(b, c);
	return dot(a, p);
}

bool Mesh::
	intersect_triangle(const Triangle &_triangle,
					   const Ray &_ray,
					   vec3 &_intersection_point,
					   vec3 &_intersection_normal,
					   double &_intersection_t,
					   vec3 &tex) const
{
	const vec3 &p0 = vertices_[_triangle.i0].position;
	const vec3 &p1 = vertices_[_triangle.i1].position;
	const vec3 &p2 = vertices_[_triangle.i2].position;

	/** 
     * - intersect _ray with _triangle
     * Rearrange `ray.origin + t*ray.dir = a*p0 + b*p1 + (1-a-b)*p2`
     * to obtain a solvable system for a, b and t using Cramer's Rule.
     * Refer to [Cramer's Rule](https://en.wikipedia.org/wiki/Cramer%27s_rule)
     */
	vec3 col0 = p1 - p0;
	vec3 col1 = p2 - p0;
	vec3 col2 = -1 * _ray.direction;

	double denominator = determ3x3(col0, col1, col2);
	vec3 d = _ray.origin - p0;

	_intersection_t = NO_INTERSECTION;

	double t = determ3x3(col0, col1, d) / denominator;
	//if t is negative => there is no intersection point
	if (t < 0)
		return false;

	double b = determ3x3(d, col1, col2) / denominator;
	if (b < 0 || b > 1)
		return false;

	double c = determ3x3(col0, d, col2) / denominator;
	double a = 1 - b - c;

	if (c < 0 || a < 0)
		return false;

	// - store ray parameter in `_intersection_t`
	_intersection_t = t;
	// - store intersection point in `_intersection_point`
	_intersection_point = _ray(_intersection_t);
	/**
    * - store normal at intersection point in `_intersection_normal`.
    * - the normal at intersection depends on the member variable `draw_mode_`
    */
	// - if the object is flat shaded, use the triangle normal
	if (draw_mode_ == FLAT)
	{
		_intersection_normal = _triangle.normal;
		// - if the object is phong shaded, use the normal interpolation formula
	}
	else if (draw_mode_ == PHONG)
	{
		const vec3 &n0 = vertices_[_triangle.i0].normal;
		const vec3 &n1 = vertices_[_triangle.i1].normal;
		const vec3 &n2 = vertices_[_triangle.i2].normal;
		_intersection_normal = normalize(a * n0 + b * n1 + c * n2);
		// interpolate texture using barycentric coordinates
		tex = a * vertices_[_triangle.i0].textures + b * vertices_[_triangle.i1].textures + c * vertices_[_triangle.i2].textures;
	}
	// - return `true` if there is an intersection with t > 0 (in front of the viewer)
	return true;
}

//=============================================================================
