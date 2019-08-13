#include "Texture.h"
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"



 Texture::Texture(const std::string &filename)
{
    int comp;
    unsigned char* image = stbi_load(filename.c_str(), &width, &height, &comp, STBI_rgb);
    if (!image){
        throw std::runtime_error("file not loaded");
    }
    data.reset((rgb*)image);
    

}
// get a pixel value as a vec3 from  UV coordinates
vec3 Texture::getPixel(double u, double v) const
{
    int u1 = v*width;
    int v1 = u*height;

    return data.get()[width*v1+u1].color();
}