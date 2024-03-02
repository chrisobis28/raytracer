#include "recursive.h"
#include "draw.h"
#include "extra.h"
#include "intersect.h"
#include "light.h"
#include <iostream>
#define XSEC_ERROR 1e-05F

// This function is provided as-is. You do not have to implement it.
// Given a range of rays, render out all rays and average the result
glm::vec3 renderRays(RenderState& state, std::span<const Ray> rays, int rayDepth)
{
    glm::vec3 L { 0.f };
    for (const auto& ray : rays) {
        L += renderRay(state, ray, rayDepth);
    }
    return L / static_cast<float>(rays.size());
}

// This method is provided as-is. You do not have to implement it.
// Given a camera ray (or secondary ray), tests for a scene intersection, and
// dependent on the results, evaluates the following functions which you must
// implement yourself:
// - `computeLightContribution()` and its submethods
// - `renderRaySpecularComponent()`, `renderRayTransparentComponent()`, `renderRayGlossyComponent()`
glm::vec3 renderRay(RenderState& state, Ray ray, int rayDepth)
{
    // Trace the ray into the scene. If nothing was hit, return early
    HitInfo hitInfo;
    if (!state.bvh.intersect(state, ray, hitInfo)) {
        drawRay(ray, glm::vec3(1, 0, 0));
        return sampleEnvironmentMap(state, ray);
    }

    // Return value: the light along the ray
    // Given an intersection, estimate the contribution of scene lights at this intersection
    glm::vec3 Lo = computeLightContribution(state, ray, hitInfo);

    // Draw an example debug ray for the incident ray (feel free to modify this for yourself)
    drawRay(ray, glm::vec3(1.0f - ((float)rayDepth / 6.0f)));

    // Given that recursive components are enabled, and we have not exceeded maximum depth,
    // estimate the contribution along these components
    if (rayDepth < 6) {
        bool isReflective = glm::any(glm::notEqual(hitInfo.material.ks, glm::vec3(0.0f)));
        bool isTransparent = hitInfo.material.transparency != 1.f;

        // Default, specular reflections
        if (state.features.enableReflections && !state.features.extra.enableGlossyReflection && isReflective) {
            renderRaySpecularComponent(state, ray, hitInfo, Lo, rayDepth);
        }

        // Alternative, glossy reflections
        if (state.features.enableReflections && state.features.extra.enableGlossyReflection && isReflective) {
            renderRayGlossyComponent(state, ray, hitInfo, Lo, rayDepth);
        }

        // Transparency passthrough
        if (state.features.enableTransparency && isTransparent) {
            renderRayTransparentComponent(state, ray, hitInfo, Lo, rayDepth);
        }
    }

    return Lo;
}

// Computes a relative error indicator for an intersection, bounded between 0 and 1,
// proportional to the square of the cosine between normal and ray direction.
// This error indicator estimates how sharp is the ray relative to the surface.
const float FLT_ERROR = FLT_EPSILON * 8.0f;

// TODO: Standard feature
// Given an incident ray and a intersection point, generate a mirrored ray
// - Ray;     the indicent ray
// - HitInfo; hit struct for the intersection point
// - return;  a reflected ray
// This method is unit-tested, so do not change the function signature.
// if you use glm::reflect, you will not get points for this method!
Ray generateReflectionRay(Ray ray, HitInfo hitInfo)
{
    Ray reflected {
        .origin = ray.origin + ray.t * ray.direction,
        .direction = ray.direction - 2.0f * glm::dot(ray.direction, hitInfo.normal) * hitInfo.normal,
        .t = FLT_MAX
    };
    // Prevent self reflection
    reflected.origin += XSEC_ERROR * reflected.direction;
#if ENABLE_MOTION_BLUR
    reflected.time = ray.time;
#endif

    // draw hit point
    drawSphere(reflected.origin, 0.03f, glm::vec3 { 1, 0, 0 });
    // draw hit normal
    drawRay(Ray { reflected.origin, hitInfo.normal, 1 }, glm::vec3 { 0, 0, 1 });
    return reflected;
}

// TODO: Standard feature
// Given an incident ray and a intersection point, generate a passthrough ray for transparency,
// starting at the intersection point and continuing in the same direction.
// - Ray;     the indicent ray
// - HitInfo; hit struct for the intersection point
// - return;  a passthrough ray for transparency
// This method is unit-tested, so do not change the function signature.
Ray generatePassthroughRay(Ray ray, HitInfo hitInfo)
{
    Ray passthrough {
        .origin = ray.origin + (ray.t + XSEC_ERROR) * ray.direction,
        .direction = ray.direction,
        .t = FLT_MAX,
    };
#if ENABLE_MOTION_BLUR
    passthrough.time = ray.time;
#endif

    // draw hit point
    drawSphere(passthrough.origin, 0.03f, glm::vec3 { 1, 0, 0 });
    // draw passthrough
    drawRay(passthrough, glm::vec3 { 0, 1, 0 });
    return passthrough;
}

// TODO: standard feature
// Given a camera ray (or secondary ray) and an intersection, evaluates the contribution
// of a mirrored ray, recursively evaluating renderRay(..., depth + 1) along this ray,
// and adding the result times material.ks to the current intersection's hit color.
// - state;    the active scene, feature config, bvh, and sampler
// - ray;      camera ray
// - hitInfo;  intersection object
// - hitColor; current color at the current intersection, which this function modifies
// - rayDepth; current recursive ray depth
// This method is unit-tested, so do not change the function signature.
void renderRaySpecularComponent(RenderState& state, Ray ray, const HitInfo& hitInfo, glm::vec3& hitColor, int rayDepth)
{
    Ray reflected = generateReflectionRay(ray, hitInfo);
    glm::vec3 Lo = renderRay(state, reflected, rayDepth + 1);
    hitColor += Lo * hitInfo.material.ks;
}

// TODO: standard feature
// Given a camera ray (or secondary ray) and an intersection, evaluates the contribution
// of a passthrough transparent ray, recursively evaluating renderRay(..., depth + 1) along this ray,
// and correctly alpha blending the result with the current intersection's hit color
// - state;    the active scene, feature config, bvh, and sampler
// - ray;      camera ray
// - hitInfo;  intersection object
// - hitColor; current color at the current intersection, which this function modifies
// - rayDepth; current recursive ray depth
// This method is unit-tested, so do not change the function signature.
void renderRayTransparentComponent(RenderState& state, Ray ray, const HitInfo& hitInfo, glm::vec3& hitColor, int rayDepth)
{
    Ray passthrough = generatePassthroughRay(ray, hitInfo);
    glm::vec3 Lo = renderRay(state, passthrough, rayDepth + 1);
    hitColor += Lo * (1.0f - hitInfo.material.transparency);
}