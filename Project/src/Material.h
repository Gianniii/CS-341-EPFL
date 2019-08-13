//=============================================================================
//                                                                            
//   Exercise code for the lecture
//   "Introduction to Computer Graphics"
//   by Prof. Dr. Mario Botsch, Bielefeld University
//
//   Copyright (C) Computer Graphics Group, Bielefeld University.
//                                                                            
//=============================================================================

#ifndef MATERIAL_H
#define MATERIAL_H


//== INCLUDES =================================================================

#include "vec3.h"
#include "Texture.h"

//== CLASS DEFINITION =========================================================


/// \class Material Material.h
/// This class stores all material parameters.
struct Material
{
    /// ambient color
    vec3   ambient;

    /// diffuse color
    vec3   diffuse;

    /// specular color
    vec3   specular;

    /// shininess factor
    double shininess;

    /// reflectivity factor (1=perfect mirror, 0=no reflection).
    double mirror;

    /// refractive index
    double refraction;

    /// add absorbtion in the material
    vec3 absorbtion;
    //optional parameter of texture
    std::unique_ptr<Texture> texture;

};


//-----------------------------------------------------------------------------


/// read material from stream
inline std::istream& operator>>(std::istream& is, Material& m)
{
    is >> m.ambient >> m.diffuse >> m.specular >> m.shininess >> m.mirror >> m.refraction >> m.absorbtion;
    return is;
}


//=============================================================================
#endif // MATERIAL_H defined
//=============================================================================
