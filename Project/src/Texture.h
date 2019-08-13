
#include "vec3.h"
#include <stdexcept>
#include <memory>
#include <string>
#include "stb_image.h"
// using namespace cimg_library;

class Texture
{

public:
    int width;
    int height;
    // struct to read independent channel from a pixel
    struct rgb
    {
        uint8_t r, g, b;
        // get a normalized pixel value as a vec3 object from a pixel
        vec3 color() const
        {
            return vec3(r/100.0, g/100.0 , b/100.0 ); 
        }
    };
    //array storing the image
    std::unique_ptr<rgb> data;

    Texture(const std::string &texturePath);

    // get a pixel value as a vec3 from  UV coordinates
    vec3 getPixel(double u, double v) const;
};
