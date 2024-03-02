#include "render.h"
#include "texture.h"
#include <cmath>
#include <fmt/core.h>
#include <glm/geometric.hpp>
#include <glm/gtx/string_cast.hpp>
#include <shading.h>

// This function is provided as-is. You do not have to implement it (unless
// you need to for some extra feature).
// Given render state and an intersection, based on render settings, sample
// the underlying material data in the expected manner.
glm::vec3 sampleMaterialKd(RenderState& state, const HitInfo& hitInfo)
{
    if (state.features.enableTextureMapping && hitInfo.material.kdTexture) {
        if (state.features.enableBilinearTextureFiltering) {
            return sampleTextureBilinear(*hitInfo.material.kdTexture, hitInfo.texCoord);
        } else {
            return sampleTextureNearest(*hitInfo.material.kdTexture, hitInfo.texCoord);
        }
    } else {
        return hitInfo.material.kd;
    }
}

// This function is provided as-is. You do not have to implement it.
// Given a camera direction, a light direction, a relevant intersection, and a color coming in
// from the light, evaluate the scene-selected shading model, returning the reflected light towards the target.
glm::vec3 computeShading(RenderState& state, const glm::vec3& cameraDirection, const glm::vec3& lightDirection, const glm::vec3& lightColor, const HitInfo& hitInfo)
{
    // Hardcoded linear gradient. Feel free to modify this
    // Modified this to be a bit more clear in testing
    static LinearGradient gradient = {
        .components = {
            { 0.1f, glm::vec3(0, 0, 1) },
            { 0.22f, glm::vec3(1, 0, 0) },
            { 0.5f, glm::vec3(0, 0, 1) },
            { 0.78f, glm::vec3(1, 0, 0) },
            { 0.9f, glm::vec3(0, 1, 0) },
        }
    };

    if (state.features.enableShading) {
        switch (state.features.shadingModel) {
        case ShadingModel::Lambertian:
            return computeLambertianModel(state, cameraDirection, lightDirection, lightColor, hitInfo);
        case ShadingModel::Phong:
            return computePhongModel(state, cameraDirection, lightDirection, lightColor, hitInfo);
        case ShadingModel::BlinnPhong:
            return computeBlinnPhongModel(state, cameraDirection, lightDirection, lightColor, hitInfo);
        case ShadingModel::LinearGradient:
            return computeLinearGradientModel(state, cameraDirection, lightDirection, lightColor, hitInfo, gradient);
        };
    }

    return lightColor * sampleMaterialKd(state, hitInfo);
}

// Given a camera direction, a light direction, a relevant intersection, and a color coming in
// from the light, evaluate a Lambertian diffuse shading, returning the reflected light towards the target.
glm::vec3 computeLambertianModel(RenderState& state, const glm::vec3& cameraDirection, const glm::vec3& lightDirection, const glm::vec3& lightColor, const HitInfo& hitInfo)
{
    // Implemented basic diffuse shading since I wished to use it...

    glm::vec3 lamb = lightColor * sampleMaterialKd(state, hitInfo) * glm::dot(lightDirection, hitInfo.normal);
    return lamb;
}

// TODO: Standard feature
// Given a camera direction, a light direction, a relevant intersection, and a color coming in
// from the light, evaluate the Phong Model returning the reflected light towards the target.
// Note: materials do not have an ambient component, so you can ignore this.
// Note: use `sampleMaterialKd` instead of material.kd to automatically forward to texture
//       sampling if a material texture is available!
//
// - state;           the active scene, feature config, and the bvh
// - cameraDirection; exitant vector towards the camera (or secondary position)
// - lightDirection;  exitant vector towards the light
// - lightColor;      the color of light along the lightDirection vector
// - hitInfo;         hit object describing the intersection point
// - return;          the result of shading along the cameraDirection vector
//
// This method is unit-tested, so do not change the function signature.
glm::vec3 computePhongModel(RenderState& state, const glm::vec3& cameraDirection, const glm::vec3& lightDirection, const glm::vec3& lightColor, const HitInfo& hitInfo)
{

    // Calculating the reflection ray and normalizing it
    glm::vec3 reflect = glm::reflect(-lightDirection, -hitInfo.normal);
    reflect = glm::normalize(reflect);
    float dotProduct = glm::dot(cameraDirection, reflect);

    // Calculating the specular shading of the Phong Model
    glm::vec3 shading = lightColor * hitInfo.material.ks * dotProduct;

    // Returning the resulting colour
    return shading + computeLambertianModel(state, cameraDirection, lightDirection, lightColor, hitInfo);
}

// TODO: Standard feature
// Given a camera direction, a light direction, a relevant intersection, and a color coming in
// from the light, evaluate the Blinn-Phong Model returning the reflected light towards the target.
// Note: materials do not have an ambient component, so you can ignore this.
// Note: use `sampleMaterialKd` instead of material.kd to automatically forward to texture
//       sampling if a material texture is available!
//
// - state;           the active scene, feature config, and the bvh
// - cameraDirection; exitant vector towards the camera (or secondary position)
// - lightDirection;  exitant vector towards the light
// - lightColor;      the color of light along the lightDirection vector
// - hitInfo;         hit object describing the intersection point
// - return;          the result of shading along the cameraDirection vector
//
// This method is unit-tested, so do not change the function signature.
glm::vec3 computeBlinnPhongModel(RenderState& state, const glm::vec3& cameraDirection, const glm::vec3& lightDirection, const glm::vec3& lightColor, const HitInfo& hitInfo)
{
    // Calculating H
    glm::vec3 H = glm::normalize(lightDirection + cameraDirection);
    float dotProduct = glm::dot(H, hitInfo.normal);

    // Calculating the specular shading of the Blinn-Phong Model
    glm::vec3 shading = lightColor * hitInfo.material.ks * dotProduct;

    // Returning the resulting colour
    return shading + computeLambertianModel(state, cameraDirection, lightDirection, lightColor, hitInfo);
}

// TODO: Standard feature
// Given a number ti between [-1, 1], sample from the gradient's components and return the
// linearly interpolated color, for which ti lies in the interval between the t-values of two
// components, or on a boundary. If ti falls outside the gradient's smallest/largest components,
// the nearest component must be sampled.
// - ti; a number between [-1, 1]
// This method is unit-tested, so do not change the function signature.
glm::vec3 LinearGradient::sample(float ti) const
{
    // Sampling the linear gradient
    LinearGradient::components;

    // Initialising low and high Values

    glm::vec3 lowColor = glm::vec3(0.0f);
    float lowT = 0.0f;

    glm::vec3 highColor = glm::vec3(1.0f);
    float highT = 1.0f;

    for (LinearGradient::Component vecj : components) {
        // Searching for closest lower bound
        if (vecj.t < ti) {
            if (vecj.t > lowT) {
                lowColor = vecj.color;
                lowT = vecj.t;
            }
        }
        // Searching for closest upper bound
        if (vecj.t > ti) {
            if (vecj.t < highT) {
                highColor = vecj.color;
                highT = vecj.t;
            }
        }
    }
    // Multiplying each factor by it's weight to return 1
    lowColor *= abs(ti - lowT) * lowColor;
    highColor *= abs(highT - ti) * highColor;

    // Returning the resulting colour
    glm::vec3 resultColor = glm::vec3 { lowColor.x + highColor.x, lowColor.y + highColor.y, lowColor.z + highColor.z };
    return resultColor;
}

// TODO: Standard feature
// Given a camera direction, a light direction, a relevant intersection, and a color coming in
// from the light, evaluate a diffuse shading model, such that the diffuse component is sampled not
// from the intersected material, but a provided linear gradient, based on the cosine of theta
// as defined in the diffuse shading part of the Phong model.
//
// - state;           the active scene, feature config, and the bvh
// - cameraDirection; exitant vector towards the camera (or secondary position)
// - lightDirection;  exitant vector towards the light
// - lightColor;      the color of light along the lightDirection vector
// - hitInfo;         hit object describing the intersection point
// - gradient;        the linear gradient object
// - return;          the result of shading
//
// This method is unit-tested, so do not change the function signature.
glm::vec3 computeLinearGradientModel(RenderState& state, const glm::vec3& cameraDirection, const glm::vec3& lightDirection, const glm::vec3& lightColor, const HitInfo& hitInfo, const LinearGradient& gradient)
{
    float cos_theta = glm::dot(lightDirection, hitInfo.normal);
    glm::vec3 kd = gradient.sample(cos_theta);
    glm::vec3 color = lightColor * kd * cos_theta;
    return color;
}