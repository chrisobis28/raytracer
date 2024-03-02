#include "light.h"
#include "bvh_interface.h"
#include "config.h"
#include "draw.h"
#include "intersect.h"
#include "render.h"
#include "scene.h"
#include "shading.h"
// Suppress warnings in third-party code.
#include <framework/disable_all_warnings.h>
DISABLE_WARNINGS_PUSH()
#include <glm/geometric.hpp>
DISABLE_WARNINGS_POP()

// TODO: Standard feature
// Given a single segment light, transform a uniformly distributed 1d sample in [0, 1),
// into a uniformly sampled position and an interpolated color on the segment light,
// and write these into the reference return values.
// - sample;    a uniformly distributed 1d sample in [0, 1)
// - light;     the SegmentLight object, see `common.h`
// - position;  reference return value of the sampled position on the light
// - color;     reference return value of the color emitted by the light at the sampled position
// This method is unit-tested, so do not change the function signature.
void sampleSegmentLight(const float& sample, const SegmentLight& light, glm::vec3& position, glm::vec3& color)
{
    // TODO: implement this function.
    // Finding the position and colour by multiplying the endpoints positions and colours by weights decided by the sample
    position = light.endpoint0 + sample * (light.endpoint1 - light.endpoint0);
    color = light.color0 + sample * (light.color1 - light.color0);
}

// TODO: Standard feature
// Given a single paralellogram light, transform a uniformly distributed 2d sample in [0, 1),
// into a uniformly sampled position and interpolated color on the paralellogram light,
// and write these into the reference return values.
// - sample;   a uniformly distributed 2d sample in [0, 1)
// - light;    the ParallelogramLight object, see `common.h`
// - position; reference return value of the sampled position on the light
// - color;    reference return value of the color emitted by the light at the sampled position
// This method is unit-tested, so do not change the function signature.
void sampleParallelogramLight(const glm::vec2& sample, const ParallelogramLight& light, glm::vec3& position, glm::vec3& color)
{
    // Finding the position and colour by multiplying the endpoints positions and colours by weights decided by the sample
    // position = light.v0 + (light.edge01 * sample.x) + light.v0 + (light.edge02 * sample.y);
    // color = (sample.x) * light.color0 / 2.0f + abs(sample.x - 1.0f) * light.color1 / 2.0f
    //    + (sample.y) * light.color2 / 2.0f + abs(sample.y - 1.0f) * light.color3 / 2.0f;
    // TODO: implement this function.
    position = light.v0 + sample.x * light.edge01 + sample.y * light.edge02;

    float xWeight = 1.0f - sample.x;
    float yWeight = 1.0f - sample.y;

    color = (xWeight * yWeight) * light.color0 + (sample.x * yWeight) * light.color1 + (xWeight * sample.y) * light.color2 + (sample.x * sample.y) * light.color3;
}

// TODO: Standard feature
// Given a sampled position on some light, and the emitted color at this position, return whether
// or not the light is visible from the provided ray/intersection.
// For a description of the method's arguments, refer to 'light.cpp'
// - state;         the active scene, feature config, and the bvh
// - lightPosition; the sampled position on some light source
// - lightColor;    the sampled color emitted at lightPosition
// - ray;           the incident ray to the current intersection
// - hitInfo;       information about the current intersection
// - return;        whether the light is visible (true) or not (false)
// This method is unit-tested, so do not change the function signature.
bool visibilityOfLightSampleBinary(RenderState& state, const glm::vec3& lightPosition, const glm::vec3& lightColor, const Ray& ray, const HitInfo& hitInfo)
{
    if (!state.features.enableShadows) {
        // Shadows are disabled in the renderer
        return true;
    } else {
        // Shadows are enabled in the renderer
        // TODO: implement this function; currently, the light simply passes through

        glm::vec3 posIntersect = ray.origin + ray.direction * ray.t;

        glm::vec3 lighttoPointRay = lightPosition - posIntersect;

        HitInfo tempInfo = HitInfo(hitInfo);
        Ray tempRay = Ray(ray);

        tempRay.direction = glm::normalize(lighttoPointRay);

        tempRay.origin = posIntersect + 0.00001f * tempRay.direction;

        state.bvh.intersect(state, tempRay, tempInfo);

        if (tempRay.t <= glm::length(lighttoPointRay)) {
            return false;
        }

        return true;
    }
}

// TODO: Standard feature
// Given a sampled position on some light, and the emitted color at this position, return the actual
// light that is visible from the provided ray/intersection, or 0 if this is not the case.
// Use the following blending operation: lightColor = lightColor * kd * (1 - alpha)
// Please reflect within 50 words in your report on why this is incorrect, and illustrate
// two examples of what is incorrect.
//
// - state;         the active scene, feature config, and the bvh
// - lightPosition; the sampled position on some light source
// - lightColor;    the sampled color emitted at lightPosition
// - ray;           the incident ray to the current intersection
// - hitInfo;       information about the current intersection
// - return;        the visible light color that reaches the intersection
//
// This method is unit-tested, so do not change the function signature.
glm::vec3 visibilityOfLightSampleTransparency(RenderState& state, const glm::vec3& lightPosition, const glm::vec3& lightColor, const Ray& ray, const HitInfo& hitInfo)
{
     glm::vec3 cameraDirection = -ray.direction;

     glm::vec3 posIntersect = ray.origin + ray.direction * ray.t;

     glm::vec3 lightToIntersection = lightPosition - posIntersect;

     if (!state.features.enableTransparency) {
         return computeShading(state, cameraDirection, lightToIntersection, lightColor, hitInfo);
     }

     HitInfo tempInfo = HitInfo(hitInfo);
     Ray tempRay = Ray(ray);

     tempRay.direction = glm::normalize(lightToIntersection);

     tempRay.origin = posIntersect + 0.00001f * tempRay.direction;

     bool intersection = state.bvh.intersect(state, tempRay, tempInfo);

     if (intersection && tempRay.t <= glm::length(lightToIntersection)) {
         return computeShading(state, cameraDirection, lightToIntersection, lightColor, hitInfo) * sampleMaterialKd(state, tempInfo) * (1.0f - tempInfo.material.transparency);
     }
     if (tempRay.t > glm::length(lightToIntersection)) {
         return glm::vec3(0.0f);
     }

    // Using the following blending operation: lightColor = lightColor * kd * (1 - alpha)
     return computeShading(state, cameraDirection, lightToIntersection, lightColor, hitInfo);

}

// TODO: Standard feature
// Given a single point light, compute its contribution towards an incident ray at an intersection point.
//
// Hint: you should use `visibilityOfLightSample()` to account for shadows, and if the light is visible, use
//       the result of `computeShading()`, whose submethods you should probably implement first in `shading.cpp`.
//
// - state;   the active scene, feature config, bvh, and a thread-safe sampler
// - light;   the PointLight object, see `common.h`
// - ray;     the incident ray to the current intersection
// - hitInfo; information about the current intersection
// - return;  reflected light along the incident ray, based on `computeShading()`
//
// This method is unit-tested, so do not change the function signature.
glm::vec3 computeContributionPointLight(RenderState& state, const PointLight& light, const Ray& ray, const HitInfo& hitInfo)
{
    // TODO: modify this function to incorporate visibility corerctly

    if (visibilityOfLightSample(state, light.position, light.color, ray, hitInfo) != glm::vec3(0.0f)) {
        return visibilityOfLightSampleTransparency(state, light.position, light.color, ray, hitInfo);
    } else if (visibilityOfLightSample(state, light.position, light.color, ray, hitInfo) == glm::vec3(0.0f) and !state.features.enableTransparency) {
        return glm::vec3(0.0f);
    } else {
        glm::vec3 p = ray.origin + ray.t * ray.direction;
        glm::vec3 l = glm::normalize(light.position - p);
        glm::vec3 v = -ray.direction;
        return computeShading(state, v, l, light.color, hitInfo);
    }
}

// TODO: Standard feature
// Given a single segment light, compute its contribution towards an incident ray at an intersection point
// by integrating over the segment, taking `numSamples` samples from the light source.
//
// Hint: you can sample the light by using `sampleSegmentLight(state.sampler.next_1d(), ...);`, which
//       you should implement first.
// Hint: you should use `visibilityOfLightSample()` to account for shadows, and if the sample is visible, use
//       the result of `computeShading()`, whose submethods you should probably implement first in `shading.cpp`.
//
// - state;      the active scene, feature config, bvh, and a thread-safe sampler
// - light;      the SegmentLight object, see `common.h`
// - ray;        the incident ray to the current intersection
// - hitInfo;    information about the current intersection
// - numSamples; the number of samples you need to take
// - return;     accumulated light along the incident ray, based on `computeShading()`
//
// This method is unit-tested, so do not change the function signature.
glm::vec3 computeContributionSegmentLight(RenderState& state, const SegmentLight& light, const Ray& ray, const HitInfo& hitInfo, uint32_t numSamples)
{
    // TODO: implement this function; repeat numSamples times:
    // - sample the segment light
    // - test the sample's visibility
    // - then evaluate the phong model
    glm::vec3 v = -ray.direction;

    glm::vec3 totalPos = glm::vec3(0.0f);
    glm::vec3 totalColor = glm::vec3(0.0f);

    HitInfo tempInfo = HitInfo(hitInfo);
    Ray tempRay = Ray(ray);

    bool intersection = state.bvh.intersect(state, tempRay, tempInfo);

    SegmentLight lightSample = light;
    glm::vec3 colori = glm::vec3(1.0f);
    glm::vec3 p = ray.origin + ray.t * ray.direction;
    for (int i = 0; i < numSamples; i++) {

        sampleSegmentLight(state.sampler.next_1d(), light, p, colori);

        PointLight lightPoint = PointLight(p, colori);

        totalColor += computeContributionPointLight(state, lightPoint, tempRay, hitInfo);
        totalPos += p;
    }
    glm::vec3 shadedColor = totalColor * (3.0f / numSamples);
    glm::vec3 weightP = totalPos * (3.0f / numSamples);

    glm::vec3 lightDir = glm::normalize(weightP);
    return computeShading(state, v, lightDir, shadedColor, hitInfo);
}

// TODO: Standard feature
// Given a single parralelogram light, compute its contribution towards an incident ray at an intersection point
// by integrating over the parralelogram, taking `numSamples` samples from the light source, and applying
// shading.
//
// Hint: you can sample the light by using `sampleParallelogramLight(state.sampler.next_1d(), ...);`, which
//       you should implement first.
// Hint: you should use `visibilityOfLightSample()` to account for shadows, and if the sample is visible, use
//       the result of `computeShading()`, whose submethods you should probably implement first in `shading.cpp`.
//
// - state;      the active scene, feature config, bvh, and a thread-safe sampler
// - light;      the ParallelogramLight object, see `common.h`
// - ray;        the incident ray to the current intersection
// - hitInfo;    information about the current intersection
// - numSamples; the number of samples you need to take
// - return;     accumulated light along the incident ray, based on `computeShading()`
//
// This method is unit-tested, so do not change the function signature.
glm::vec3 computeContributionParallelogramLight(RenderState& state, const ParallelogramLight& light, const Ray& ray, const HitInfo& hitInfo, uint32_t numSamples)
{
    // TODO: implement this function; repeat numSamples times:
    // - sample the parallellogram light
    // - test the sample's visibility
    // - then evaluate the phong model
    glm::vec3 v = -ray.direction;

    glm::vec3 totalPos = glm::vec3(0.0f);
    glm::vec3 totalColor = glm::vec3(0.0f);

    HitInfo tempInfo = HitInfo(hitInfo);
    Ray tempRay = Ray(ray);

    bool intersection = state.bvh.intersect(state, tempRay, tempInfo);

    ParallelogramLight lightSample = light;
    glm::vec3 colori = glm::vec3(1.0f);
    glm::vec3 p = ray.origin + ray.t * ray.direction;
    for (int i = 0; i < numSamples; i++) {

        sampleParallelogramLight(state.sampler.next_2d(), light, p, colori);

        PointLight lightPoint = PointLight(p, colori);

        totalColor += computeContributionPointLight(state, lightPoint, tempRay, hitInfo);
        totalPos += p;
    }
    glm::vec3 shadedColor = totalColor * (3.0f / numSamples);
    glm::vec3 weightP = totalPos * (3.0f / numSamples);

    glm::vec3 lightDir = glm::normalize(weightP);
    return computeShading(state, v, lightDir, shadedColor, hitInfo);
}

// This function is provided as-is. You do not have to implement it.
// Given a sampled position on some light, and the emitted color at this position, return the actual
// light that is visible from the provided ray/intersection, or 0 if this is not the case.
// This forowards to `visibilityOfLightSampleBinary`/`visibilityOfLightSampleTransparency` based on settings.
//
// - state;         the active scene, feature config, and the bvh
// - lightPosition; the sampled position on some light source
// - lightColor;    the sampled color emitted at lightPosition
// - ray;           the incident ray to the current intersection
// - hitInfo;       information about the current intersection
// - return;        the visible light color that reaches the intersection
//
// This method is unit-tested, so do not change the function signature.
glm::vec3 visibilityOfLightSample(RenderState& state, const glm::vec3& lightPosition, const glm::vec3& lightColor, const Ray& ray, const HitInfo& hitInfo)
{
    if (!state.features.enableShadows) {
        // Shadows are disabled in the renderer
        return lightColor;
    } else if (!state.features.enableTransparency) {
        // Shadows are enabled but transparency is disabled
        return visibilityOfLightSampleBinary(state, lightPosition, lightColor, ray, hitInfo) ? lightColor : glm::vec3(0);
    } else {
        // Shadows and transparency are enabled
        return visibilityOfLightSampleTransparency(state, lightPosition, lightColor, ray, hitInfo);
    }
}

// This function is provided as-is. You do not have to implement it.
glm::vec3 computeLightContribution(RenderState& state, const Ray& ray, const HitInfo& hitInfo)
{
    // Iterate over all lights
    glm::vec3 Lo { 0.0f };
    for (const auto& light : state.scene.lights) {
        if (std::holds_alternative<PointLight>(light)) {
            Lo += computeContributionPointLight(state, std::get<PointLight>(light), ray, hitInfo);
        } else if (std::holds_alternative<SegmentLight>(light)) {
            Lo += computeContributionSegmentLight(state, std::get<SegmentLight>(light), ray, hitInfo, state.features.numShadowSamples);
        } else if (std::holds_alternative<ParallelogramLight>(light)) {
            Lo += computeContributionParallelogramLight(state, std::get<ParallelogramLight>(light), ray, hitInfo, state.features.numShadowSamples);
        }
    }
    return Lo;
}