#include "texture.h"
#include "render.h"
#include <framework/image.h>

// Transform [i, j] picture coordinates into index in the image array
// Use texture extend at edges of the image
glm::vec3 getPixel(int x, int y, const Image& image) {
    x = glm::min(glm::max(0, x), image.width - 1);
    y = glm::min(glm::max(0, y), image.height - 1);
    return image.pixels[image.width * (image.height - 1 - y) + x];
}

// TODO: Standard feature
// Given an image, and relevant texture coordinates, sample the texture s.t.
// the nearest texel to the coordinates is acquired from the image.
// - image;    the image object to sample from.
// - texCoord; sample coordinates, generally in [0, 1]
// - return;   the nearest corresponding texel
// This method is unit-tested, so do not change the function signature.
glm::vec3 sampleTextureNearest(const Image& image, const glm::vec2& texCoord)
{
    // TODO: implement this function.
    // Note: the pixels are stored in a 1D array, row-major order. You can convert from (i, j) to
    //       an index using the method seen in the lecture.
    // Note: the center of the first pixel should be at coordinates (0.5, 0.5)
    // Given texcoords, return the corresponding pixel of the image
    // The pixel are stored in a 1D array of row major order
    // you can convert from position (i,j) to an index using the method seen in the lecture
    // Note, the center of the first pixel is at image coordinates (0.5, 0.5)

    int x = (int) (texCoord.x * (float) image.width);
    int y = (int) (texCoord.y * (float) image.height);
    return getPixel(x, y, image);
}

// TODO: Standard feature
// Given an image, and relevant texture coordinates, sample the texture s.t.
// a bilinearly interpolated texel is acquired from the image.
// - image;    the image object to sample from.
// - texCoord; sample coordinates, generally in [0, 1]
// - return;   the filter of the corresponding texels
// This method is unit-tested, so do not change the function signature.
glm::vec3 sampleTextureBilinear(const Image& image, const glm::vec2& texCoord)
{
    // TODO: implement this function.
    // Note: the pixels are stored in a 1D array, row-major order. You can convert from (i, j) to
    //       an index using the method seen in the lecture.
    // Note: the center of the first pixel should be at coordinates (0.5, 0.5)
    // Given texcoords, return the corresponding pixel of the image
    // The pixel are stored in a 1D array of row major order
    // you can convert from position (i,j) to an index using the method seen in the lecture
    // Note, the center of the first pixel is at image coordinates (0.5, 0.5)

    // Get the upper-left-most texel of the closest 4
    // [a]  b
    //  c   d
    int x = (int) (texCoord.x * (float) image.width - 0.5f);
    int y = (int) (texCoord.y * (float) image.height - 0.5f);

    glm::vec3 a = getPixel(x, y, image);
    glm::vec3 b = getPixel(x + 1, y, image);
    glm::vec3 c = getPixel(x, y + 1, image);
    glm::vec3 d = getPixel(x + 1, y + 1, image);

    float alpha = texCoord.x * (float) image.width - (float) x - 0.5f;
    float beta = texCoord.y * (float) image.height - (float) y - 0.5f;

    return beta * (alpha * d + (1.f - alpha) * c) + (1.f - beta) * (alpha * b + (1.f - alpha) * a);
}