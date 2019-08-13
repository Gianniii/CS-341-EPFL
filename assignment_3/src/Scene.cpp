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
#include "Scene.h"

#include "Plane.h"
#include "Sphere.h"
#include "Cylinder.h"
#include "Mesh.h"

#include <limits>
#include <map>
#include <functional>
#include <stdexcept>

#if HAS_TBB
#include <tbb/tbb.h>
#include <tbb/parallel_for.h>
#endif

//-----------------------------------------------------------------------------

Image Scene::render()
{
    // allocate new image.
    Image img(camera.width, camera.height);

    // Function rendering a full column of the image
    auto raytraceColumn = [&img, this](int x) {
        for (int y=0; y<int(camera.height); ++y)
        {
            Ray ray = camera.primary_ray(x,y);

            // compute color by tracing this ray
            vec3 color = trace(ray, 0);

            // avoid over-saturation
            color = min(color, vec3(1, 1, 1));

            // store pixel color
            img(x,y) = color;
        }
    };

    // If possible, raytrace image columns in parallel. We use TBB if available
    // and try OpenMP otherwise. Note that OpenMP only works on the latest
    // clang compilers, so macOS users will probably have the best luck with TBB.
    // You can install TBB with MacPorts/Homebrew, or from Intel:
    // https://github.com/01org/tbb/releases
#if HAS_TBB
    tbb::parallel_for(tbb::blocked_range<int>(0, camera.width), [&raytraceColumn](const tbb::blocked_range<int> &range) {
        for (size_t i = range.begin(); i < range.end(); ++i)
            raytraceColumn(i);
    });
#else
#if defined(_OPENMP)
#pragma omp parallel for
#endif
    for (int x=0; x<int(camera.width); ++x)
        raytraceColumn(x);
#endif

    // Note: compiler will elide copy.
    return img;
}

//-----------------------------------------------------------------------------

vec3 Scene::trace(const Ray& _ray, int _depth)
{
    // stop if recursion depth (=number of reflections) is too large
    if (_depth > max_depth) return vec3(0,0,0);

    // Find first intersection with an object. If an intersection is found,
    // it is stored in object, point, normal, and t.
    Object_ptr  object;
    vec3        point;
    vec3        normal;
    double      t;
    if (!intersect(_ray, object, point, normal, t))
    {
        return background;
    }

    // compute local Phong lighting (ambient+diffuse+specular)
    vec3 color = lighting(point, normal, -_ray.direction, object->material);

    const Material mat = object->material;

    /**
     * Compute reflections by recursive ray tracing:
     * - check whether `object` is reflective by checking its `material.mirror`
     * - check recursion depth
     * - generate reflected ray, compute its color contribution, and mix it with
     * the color computed by local Phong lighting (use `object->material.mirror` as weight)
     * - check whether your recursive algorithm reflects the ray `max_depth` times
     */


    if(mat.mirror > 0 && _depth <= max_depth){
        vec3 r = reflect(normalize(_ray.direction), normal); //create a vector with the direction of the reflect one
        Ray reflRay = Ray(point + r * 0.001, r); //launch a ray from the intersection point
        vec3 colorRefl = trace(reflRay, ++_depth);//doing it recursively
        return colorRefl * mat.mirror + (1 - mat.mirror) * color;//adjusting the color

    }  else {
        return color; //at the end when there is no more reflect return color
    }


}

//-----------------------------------------------------------------------------

bool Scene::intersect(const Ray& _ray, Object_ptr& _object, vec3& _point, vec3& _normal, double& _t)
{
    double  t, tmin(Object::NO_INTERSECTION);
    vec3    p, n;

    for (const auto &o: objects) // for each object
    {
        if (o->intersect(_ray, p, n, t)) // does ray intersect object?
        {
            if (t < tmin) // is intersection point the currently closest one?
            {
                tmin = t;
                _object = o.get();
                _point  = p;
                _normal = n;
                _t      = t;
            }
        }
    }

    return (tmin != Object::NO_INTERSECTION);
}

vec3 Scene::lighting(const vec3& _point, const vec3& _normal, const vec3& _view, const Material& _material)
{

     /**
     * Compute the Phong lighting:
     * - start with global ambient contribution
     * - for each light source (stored in vector `lights`) add diffuse and specular contribution
     * - only add diffuse and specular light if object is not in shadow
     */
    
    vec3 color = ambience * _material.ambient;// computing the ambient light
    
    for (Light &light : lights)
    {
        const vec3 l = normalize((light.position - _point));//creating a normal vector from the intrsection point to the light
        const vec3 n = normalize(_normal);//normalizing the normal (if it is not)
        const vec3 r = normalize(2 * n * dot(n, l) - l);//create the ray that reflect the light
        const vec3 v = _view;
        if (dot(n, l) >= 0){//checking if there is diffusion 
            vec3 specRefl = light.color * _material.diffuse * dot(n, l);//adding the diffusion a variable
            if (dot(r,v) >= 0){// checking if there is specular
                specRefl += light.color * _material.specular * pow(dot(r, v), _material.shininess);//adding the specular to the diffusion
            }
            const Ray shadowRay = Ray(_point + l * 0.001 , l);//creating a ray that go from the intersection point to the light

            double  t;
            vec3 p, normal;
            bool shadowed = false;
            for (const auto &o: objects) // for each object
            {
                if (o->intersect(shadowRay, p, normal, t)) // does this ray intersect object?
                {
                    //if yes check that the object it intersects is between the light and the object and not behind the light
                    if(norm(p - _point) <= norm(light.position - _point)){
                       shadowed = true;//if yes then it is in a shadow
                        break; 
                    }   
                }
            }
            if(!shadowed){
            // if the object is not in the shadow we add the diffuse and specular color to the ambient color
            // if not we are not adding it
                color += specRefl;
            }
        }
    }


    return color;//at the end we return the color
}

//-----------------------------------------------------------------------------

void Scene::read(const std::string &_filename)
{
    std::ifstream ifs(_filename);
    if (!ifs)
        throw std::runtime_error("Cannot open file " + _filename);

    const std::map<std::string, std::function<void(void)>> entityParser = {
        {"depth",      [&]() { ifs >> max_depth; }},
        {"camera",     [&]() { ifs >> camera; }},
        {"background", [&]() { ifs >> background; }},
        {"ambience",   [&]() { ifs >> ambience; }},
        {"light",      [&]() { lights .emplace_back(ifs); }},
        {"plane",      [&]() { objects.emplace_back(new    Plane(ifs)); }},
        {"sphere",     [&]() { objects.emplace_back(new   Sphere(ifs)); }},
        {"cylinder",   [&]() { objects.emplace_back(new Cylinder(ifs)); }},
        {"mesh",       [&]() { objects.emplace_back(new     Mesh(ifs, _filename)); }}
    };

    // parse file
    std::string token;
    while (ifs && (ifs >> token) && (!ifs.eof())) {
        if (token[0] == '#') {
            ifs.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            continue;
        }

        if (entityParser.count(token) == 0)
            throw std::runtime_error("Invalid token encountered: " + token);
        entityParser.at(token)();
    }
}


//=============================================================================
