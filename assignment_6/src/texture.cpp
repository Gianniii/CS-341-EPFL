//=============================================================================
//
//   Exercise code for the lecture "Introduction to Computer Graphics"
//     by Prof. Mario Botsch, Bielefeld University
//
//   Copyright (C) by Computer Graphics Group, Bielefeld University
//
//=============================================================================

#include "texture.h"
#include <iostream>
#include <cassert>
#include <algorithm>
#include "lodepng.h"
#include <math.h>  

//=============================================================================


Texture::Texture() :
    id_(0)
{
}


//-----------------------------------------------------------------------------


Texture::~Texture()
{
    if (id_) glDeleteTextures(1, &id_);
}


//-----------------------------------------------------------------------------

bool Texture::loadPNG(const char* filename)
{
    std::cout << "Load texture " << filename << "\n" << std::flush;

    std::vector<unsigned char> img;
    unsigned width, height;

    unsigned error = lodepng::decode(img, width, height, filename);
    if (error) {
        std::cout << "read error: " << lodepng_error_text(error) << std::endl;
        return false;
    }

    return uploadImage(img, width, height);
}

//-----------------------------------------------------------------------------


bool Texture::uploadImage(std::vector<unsigned char> &img, unsigned width, unsigned height)
{
    if (!id_) {
        std::cerr << "Texture: initialize before loading!\n";
        return false;
    }

    // flip vertically in order to adhere to how OpenGL interpretes image data
    if (height % 2) throw std::runtime_error("Image must be of even height");
    assert(height % 2 == 0);
    for (unsigned int y = 0; y < height/2; ++y) {
        for (unsigned int x = 0; x < width; ++x) {
            for (unsigned int c = 0; c < 4; ++c) {
                std::swap(img[(              y  * width + x) * 4 + c],
                          img[((height - y - 1) * width + x) * 4 + c]);
            }
        }
    }

    // upload texture data
    glActiveTexture(unit_);
    glBindTexture(type_, id_);
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
    glPixelStorei(GL_PACK_ALIGNMENT, 1);
    glTexImage2D(type_, 0, GL_RGBA, width, height, 0, GL_RGBA, GL_UNSIGNED_BYTE, &img[0]);


    if(minfilter_==GL_LINEAR_MIPMAP_LINEAR)
    {
        // comment out to disable mipmaps
        glGenerateMipmap(GL_TEXTURE_2D);
    }

    return true;
}


//-----------------------------------------------------------------------------

bool Texture::createSunBillboardTexture()
{
    std::cout << "creating sun billboard " << "\n" << std::flush;

    std::vector<unsigned char> img;
    int width = 900;
    int height = 900;
    img.resize(width*height * 4);
    float middle_width = width/2.0;
    float middle_height = height/2.0;
    float square_radius = 150 * 150;
    float square_max_opacity = 300 * 300;

    /**  Set up the texture for the sun billboard.
    *   - Draw an opaque circle with a 150 pixel radius in its middle
    *   - Outside that circle the texture should become more and more transparent to mimic a nice glow effect
    *   - Make sure that your texture is fully transparent at its borders to avoid seeing visible edges
    *   - Experiment with the color and with how fast you change the transparency until the effect satisfies you
    **/

    for (int col = 0; col < width; ++col) {
        for (int row = 0; row < height; ++row) {
            float square_distance_center = (col - middle_width) * (col - middle_width) + (row - middle_height) * (row - middle_height);
            if (square_distance_center <= square_radius){
                img[(row * width + col) * 4 + 3] = 255; // A
            } else if (square_distance_center >= square_max_opacity){
                img[(row * width + col) * 4 + 3] = 0; // A
            }else{
                img[(row * width + col) * 4 + 3] = int(255 * (1 - (square_distance_center - square_radius)/(square_max_opacity - square_radius))); // A
            }

            float dist = sqrt(pow((row - height/2),2) + pow((col - width/2),2));
            
            img[(row * width + col) * 4 + 0] = 253; // R
            img[(row * width + col) * 4 + 1] = 106; // G
            img[(row * width + col) * 4 + 2] = 2; // B
            img[(row * width + col) * 4 + 3] = dist < 200 ? (dist < 150 ? 255 : (
                150*(1 - (dist - 150)*3/150))) : 0; // A a basic linear from 150 px to 200 with transparency that goes from 255 to 0 
        }
    }

    return uploadImage(img, width, height);
}


//-----------------------------------------------------------------------------


void Texture::init(GLenum unit, GLenum type, GLint minfilter, GLint magfilter, GLint wrap)
{
    // remember this
    unit_ = unit;
    type_ = type;
    minfilter_ = minfilter;

    // activate texture unit
    glActiveTexture(unit_);

    // create texture object
    glGenTextures(1, &id_);
    glBindTexture(type_, id_);

    // set texture parameters
    glTexParameteri(type_, GL_TEXTURE_MAG_FILTER, magfilter);
    glTexParameteri(type_, GL_TEXTURE_MIN_FILTER, minfilter);
    glTexParameteri(type_, GL_TEXTURE_WRAP_S, wrap);
    glTexParameteri(type_, GL_TEXTURE_WRAP_T, wrap);
}


//-----------------------------------------------------------------------------


void Texture::bind()
{
    assert(id_);
    glActiveTexture(unit_);
    glBindTexture(type_, id_);
}



//=============================================================================
