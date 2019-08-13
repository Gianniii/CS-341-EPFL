#include "Image.h"
#include <iostream>
#include <algorithm>
#include <lodepng.h>

static bool check_ext(const std::string &path, const std::string &ext)
{
    if (ext.size() > path.size()) {
        return false;
    }
    // reading from the end
    return std::equal(ext.rbegin(), ext.rend(), path.rbegin());
}

bool Image::write(const std::string &_filename) const {
    if (check_ext(_filename, ".png")) return write_png(_filename);
    if (check_ext(_filename, ".tga")) return write_tga(_filename);

    std::cerr << "No encoder for file name " << _filename << std::endl;
    return false;
}

bool Image::write_tga(const std::string &_filename) const
{
    std::ofstream file(_filename, std::fstream::binary);
    if (!file) return false;

    file.put(0); //id length
    file.put(0); //no color map
    file.put(2); //uncompressed image
    file.put(0); //offset color map table
    file.put(0); //
    file.put(0); //number of entries
    file.put(0); //
    file.put(0); //bits per pixel
    file.put(0); //abs coordinate lower left display in x direction
    file.put(0); //
    file.put(0); //abs coordinate lower left display in y direction
    file.put(0); //
    file.put((width_  & 0x00FF)); //width in pixels
    file.put((width_  & 0xFF00)/ 256);
    file.put((height_ & 0x00FF)); //height in pixels
    file.put((height_ & 0xFF00)/ 256);
    file.put(24); //bits per pixel
    file.put(0); //image descriptor

    for (const vec3 &color : pixels_)
    {
        file.put(static_cast<unsigned char>(255.0 * color[2]));
        file.put(static_cast<unsigned char>(255.0 * color[1]));
        file.put(static_cast<unsigned char>(255.0 * color[0]));
    }

    file.close();
    return true;
}

bool Image::write_png(const std::string &_filename) const {

    std::vector<uint8_t> image_data(pixels_.size() * 3);

    for (unsigned int y = 0; y < height(); ++y) {
        for (unsigned int x = 0; x < width(); ++x) {
            unsigned int row = height() - 1 - y; // flip vertically (PNG origin is upper left)
            for (uint8_t c = 0; c < 3; ++c)
                image_data[3 * (row * width() + x) + c] = static_cast<unsigned char>(255.0 * (*this)(x, y)[c]);
        }
    }

    return lodepng::encode(_filename, image_data, width(), height(), LCT_RGB);
}
