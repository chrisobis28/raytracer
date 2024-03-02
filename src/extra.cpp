#include "extra.h"
#include "bvh.h"
#include "light.h"
#include "recursive.h"
#include "shading.h"
#include "draw.h"
#include <framework/trackball.h>
#include "texture.h"
#include "intersect.h"
#include "draw.h"

// TODO; Extra feature
// Given the same input as for `renderImage()`, instead render an image with your own implementation
// of Depth of Field. Here, you generate camera rays s.t. a focus point and a thin lens camera model
// are in play, allowing objects to be in and out of focus.
// This method is not unit-tested, but we do expect to find it **exactly here**, and we'd rather
// not go on a hunting expedition for your implementation, so please keep it here!
void renderImageWithDepthOfField(const Scene& scene, const BVHInterface& bvh, const Features& features, const Trackball& camera, Screen& screen)
{
    if (!features.extra.enableDepthOfField) {
        return;
    }

    // ...
}

// *
// *
// MOTION BLUR FUNCTIONS
// *
// Disabled by default, change the value of flag "ENABLE_MOTION_BLUR" in mesh.h to include motion blur
// source code in the compilation.
// *
// *
#if ENABLE_MOTION_BLUR

// Returns displacement of the curve given time t
// The curve is defined by Bézier curve with p0 in (0, 0, 0)
// p1, p2 and p3 are passed to the function, p1 and p2 have to be in a cube bounded by p0 and p3
glm::vec3 cubicBezier(glm::vec3 p1, glm::vec3 p2, glm::vec3 p3, float t) {
    float t_inv = 1.f - t;
    return 3.f * t_inv * t_inv * t * p1 + 3.f * t_inv * t * t * p2 + t * t * t * p3;
}

glm::vec3 quadBezier(glm::vec3 p1, glm::vec3 p2, float t) {
    return 2.f * (1.f - t) * t * p1 + t * t * p2;
}

// Returns new vector that is a copy of the original, transformed over time
glm::vec3 applyTransform(glm::vec3 point, const Transform& transform, float time)
{
    switch (transform.type) {
    case LINEAR: {
        return point + time * transform.data[0];
    }
    case BEZIER: {
        return point + cubicBezier(transform.data[0], transform.data[1], transform.data[2], time);
    }
    case ROTATION: {
        // TODO: Implement rotation
        return point;
    }
    default: {
        return point;
    }
    }
}

// Returns new vertex that is a copy of the original, transformed over time
Vertex applyTransform(const Vertex& v, const Transform& transform, float time) {
    Vertex vTransformed = v;
    vTransformed.position = applyTransform(vTransformed.position, transform, time);
    return vTransformed;
}

// Adds transform functions to scene meshes dependent on the scene type
void updateSceneWithTransforms(Scene& scene) {
    switch (scene.type) {
    case Lambo: {
        Transform bezier {
            BEZIER,
            { glm::vec3 {-7, 0, 10}, glm::vec3 {-5, 0, 20}, glm::vec3 {0, 0, 30} }
        };
        scene.meshes[0].transform = bezier;
        scene.meshes[1].transform = bezier;
    } break;
    case Spheres: {
        scene.spheres[0].transform = Transform {
            BEZIER,
            { glm::vec3 { 0, 2, 1 }, glm::vec3 { 3, 1, 2 }, glm::vec3 { 3, 3, 3 } }
        };
        scene.spheres[1].transform = Transform {
            BEZIER,
            { glm::vec3 { 0, -2, -1 }, glm::vec3 { -3, -1, -2 }, glm::vec3 { -3, -3, -3 } }
        };
        scene.spheres[2].transform = Transform {
            BEZIER,
            { glm::vec3 { 0, 6, 3 }, glm::vec3 { 9, 3, 6 }, glm::vec3 { 9, 9, 9 } }
        };
    } break;
    default: {
        Transform bezier {
            BEZIER,
            { glm::vec3 { 0, 2, 1 }, glm::vec3 { 3, 1, 2 }, glm::vec3 { 3, 3, 3 } }
        };
        for (auto& mesh : scene.meshes) {
            mesh.transform = bezier;
        }
        for (auto& sphere : scene.spheres) {
            sphere.transform = bezier;
        }
    } break;
    }
}

void drawTransform(const Transform& transform, const glm::vec3& p0, float time) {
    std::vector<glm::vec3> points { p0 };
    points.resize(100);
    float time_step = time * 0.01f;
    for (int i = 1; i < 100; i++) {
        points[i] = applyTransform(p0, transform, (float) i * time_step);
    }
    drawLineSegment(points);
}
#endif

// TODO; Extra feature
// Given the same input as for `renderImage()`, instead render an image with your own implementation
// of motion blur. Here, you integrate over a time domain, and not just the pixel's image domain,
// to give objects the appearance of "fast movement".
// This method is not unit-tested, but we do expect to find it **exactly here**, and we'd rather
// not go on a hunting expedition for your implementation, so please keep it here!
void renderImageWithMotionBlur(const Scene& scene, const BVHInterface& bvh, const Features& features, const Trackball& camera, Screen& screen)
{
    if (!features.extra.enableMotionBlur) return;

#if ENABLE_MOTION_BLUR
    uint32_t samples = features.extra.numTimeSamples;
    float shutter_speed = features.extra.shutterSpeed;

#ifdef NDEBUG // Enable multi threading in Release mode
#pragma omp parallel for schedule(guided)
#endif
    for (int y = 0; y < screen.resolution().y; y++) {
        for (int x = 0; x != screen.resolution().x; x++) {
            // Assemble useful objects on a per-pixel basis; e.g. a per-thread sampler
            // Note; we seed the sampler for consistent behavior across frames
            RenderState state = {
                .scene = scene,
                .features = features,
                .bvh = bvh,
                .sampler = { static_cast<uint32_t>(screen.resolution().y * x + y) }
            };

            // Sample color at random timestamps in the shutter speed interval
            glm::vec3 L(0);
            for (uint32_t i = 0; i < samples; i++) {
                float timeSample = state.sampler.next_1d() * shutter_speed;
                auto rays = generatePixelRays(state, camera, { x, y }, screen.resolution());
                // Update rays time
                for (auto &ray : rays) {
                    ray.time = timeSample;
                }
                L += renderRays(state, rays);
            }
            screen.setPixel(x, y, L / (float) samples);
        }
    }

#endif
}

int factorial(const int base, const int exponent)
{
    int newExp = exponent - 1;
    if (exponent == 0 || exponent == 1) {
        return base;
    } else {
        return factorial(base * exponent, newExp);
    }
}

// TODO; Extra feature
// Given a rendered image, compute and apply a bloom post-processing effect to increase bright areas.
// This method is not unit-tested, but we do expect to find it **exactly here**, and we'd rather
// not go on a hunting expedition for your implementation, so please keep it here!
void postprocessImageWithBloom(const Scene& scene, const Features& features, const Trackball& camera, Screen& image)
{
    if (!features.extra.enableBloomEffect) {
        return;
    }

    int width = image.resolution().x;
    int height = image.resolution().y;

    float totalValue;
    int size = image.pixels().size();

    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            int index = image.indexAt(j, i);
            glm::vec3 color = image.pixels()[index] * float(factorial(index, index)) * (1.0f / float(factorial(size, size)) * float(factorial(index - size, index - size)));
            image.setPixel(j, i, color);
            totalValue = totalValue + glm::length(color);
        }
    }

    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            int index = image.indexAt(i, j);
            glm::vec3 color = image.pixels()[index] * float(1.0f / totalValue);
            image.setPixel(j, i, color);
        }
    }

    for (int i = 0; i < width; i++) {
        for (int j = 0; j < height; j++) {
            int index = image.indexAt(j, i);
            glm::vec3 color = image.pixels()[index] * float(factorial(index, index)) * (1.0f / float(factorial(size, size)) * float(factorial(index - size, index - size)));
            image.setPixel(j, i, color);
            totalValue = totalValue + glm::length(color);
        }
    }

    for (int i = 0; i < width; i++) {
        for (int j = 0; j < height; j++) {
            int index = image.indexAt(i, j);
            glm::vec3 color = image.pixels()[index] * float(1.0f / totalValue);
            image.setPixel(j, i, color);
        }
    }

    // ...
}


// TODO; Extra feature
// Given a camera ray (or reflected camera ray) and an intersection, evaluates the contribution of a set of
// glossy reflective rays, recursively evaluating renderRay(..., depth + 1) along each ray, and adding the
// results times material.ks to the current intersection's hit color.
// - state;    the active scene, feature config, bvh, and sampler
// - ray;      camera ray
// - hitInfo;  intersection object
// - hitColor; current color at the current intersection, which this function modifies
// - rayDepth; current recursive ray depth
// This method is not unit-tested, but we do expect to find it **exactly here**, and we'd rather
// not go on a hunting expedition for your implementation, so please keep it here!
void renderRayGlossyComponent(RenderState& state, Ray ray, const HitInfo& hitInfo, glm::vec3& hitColor, int rayDepth)
{
    // Generate an initial specular ray, and base secondary glossies on this ray
    // auto numSamples = state.features.extra.numGlossySamples;
    // ...
    float limit = (float)state.features.extra.numGlossySamples;

    // compute secondary rays
    Ray reflected = generateReflectionRay(ray, hitInfo);

    glm::vec3 u = glm::normalize(glm::vec3(0, 1, -reflected.direction.y / reflected.direction.z));
    glm::vec3 v = glm::normalize(glm::cross(reflected.direction, u));
    glm::vec3 intersection = ray.origin + ray.t * ray.direction;

    glm::vec3 sum = {0.0f, 0.0f, 0.0f};

    for (int i = 0; i < state.features.extra.numGlossySamples; ++i) {
        glm::vec2 x = state.sampler.next_2d();

        float polar_angle = glm::pi<float>() * 2.0f * x.x;
        float r = hitInfo.material.shininess / 64.0f * (glm::sqrt(x.x * x.x + x.y * x.y));

        float indexu = r * glm::cos(polar_angle);
        float indexv = r * glm::sin(polar_angle);

        Ray glossySample;
        glossySample.direction = glm::normalize(reflected.origin + reflected.direction + indexu * u + indexv * v);
        glossySample.origin = intersection + 0.001f * glossySample.direction;
        if (glm::dot(glossySample.direction, hitInfo.normal) > 0) {
            glm::vec3 color = renderRay(state, glossySample, rayDepth + 1);
            if (state.features.extra.enableGlossyDebug)
                drawRay(glossySample, color);
            sum += color * hitInfo.material.ks;
        }
    }

    hitColor += sum / limit;
}

// TODO; Extra feature
// Given a camera ray (or reflected camera ray) that does not intersect the scene, evaluates the contribution
// along the ray, originating from an environment map. You will have to add support for environment textures
// to the Scene object, and provide a scene with the right data to supply this.
// - state; the active scene, feature config, bvh, and sampler
// - ray;   ray object
// This method is not unit-tested, but we do expect to find it **exactly here**, and we'd rather
// not go on a hunting expedition for your implementation, so please keep it here!
glm::vec3 sampleEnvironmentMap(RenderState& state, Ray ray)
{
    if (state.features.extra.enableEnvironmentMap) {

        AxisAlignedBox cube = AxisAlignedBox(glm::vec3(-1.0f, -1.0f, -1.0f), glm::vec3(1.0f, 1.0f, 1.0f));
        Ray new_ray;
        new_ray.direction = glm::normalize(ray.direction);
        bool intersection_boolean = intersectRayWithShape(cube, new_ray);
        if (intersection_boolean) {
            glm::vec3 point = new_ray.origin + new_ray.t * new_ray.direction;

            float absX = glm::abs(point.x);
            float absY = glm::abs(point.y);
            float absZ = glm::abs(point.z);
            int faceIndex = 0;
            float v = 0, u = 0;

            if (point.x > 0 && absX >= absY && absX >= absZ) {
                faceIndex = 0;
                v = (point.y / point.x + 1) * 0.5f;
                u = (-point.z / point.x + 1) * 0.5f;
            } else if (point.x < 0 && absX >= absY && absX >= absZ) {
                faceIndex = 1;
                v = (point.y / glm::abs(point.x) + 1) * 0.5f;
                u = (point.z / glm::abs(point.x) + 1) * 0.5f;
            } else if (point.y > 0 && absY >= absX && absY >= absZ) {
                faceIndex = 2;
                v = (-point.z / point.y + 1) * 0.5f;
                u = (point.x / point.y + 1) * 0.5f;
            } else if (point.y < 0 && absY >= absX && absY >= absZ) {
                faceIndex = 3;
                v = (point.z / glm::abs(point.y) + 1) * 0.5f;
                u = (point.x / glm::abs(point.y) + 1) * 0.5f;
            } else if (point.z > 0 && absZ >= absX && absZ >= absY) {
                faceIndex = 4;
                v = (point.y / point.z + 1) * 0.5f;
                u = (point.x / point.z + 1) * 0.5f;
            } else if (point.z < 0 && absZ >= absX && absZ >= absY) {
                faceIndex = 5;
                v = (point.y / glm::abs(point.z) + 1) * 0.5f;
                u = (-point.x / glm::abs(point.z) + 1) * 0.5f;
            }

            return sampleTextureNearest(*state.scene.images[faceIndex], glm::vec2(u, v));
        } else
            return glm::vec3(0.0f);
    } else
        return glm::vec3(0.0f);
}


// TODO: Extra feature
// As an alternative to `splitPrimitivesByMedian`, use a SAH+binning splitting criterion. Refer to
// the `Data Structures` lecture for details on this metric.
// - aabb;       the axis-aligned bounding box around the given triangle set
// - axis;       0, 1, or 2, determining on which axis (x, y, or z) the split must happen
// - primitives; the modifiable range of triangles that requires splitting
// - return;     the split position of the modified range of triangles
// This method is unit-tested, so do not change the function signature.

float surfaceArea(AxisAlignedBox aabb)
{
    if (aabb.lower.x == FLT_MAX)
        return 0;
    return 2 * ((aabb.upper.x - aabb.lower.x) * (aabb.upper.y - aabb.lower.y) + (aabb.upper.x - aabb.lower.x) * (aabb.upper.z - aabb.lower.z) + (aabb.upper.z - aabb.lower.z) * (aabb.upper.y - aabb.lower.y));
}

size_t splitPrimitivesBySAHBin(const AxisAlignedBox& aabb, uint32_t axis, std::span<BVH::Primitive> primitives)
{
    using Primitive = BVH::Primitive;

    float axis_length = aabb.upper[axis] - aabb.lower[axis];
    int k = 16;

    struct Bucket {
        size_t count = 0;
        AxisAlignedBox aabb = AxisAlignedBox({ FLT_MAX, FLT_MAX, FLT_MAX }, { -FLT_MAX, -FLT_MAX, -FLT_MAX });
    };

    std::vector<Bucket> buckets(k);
    std::vector<AxisAlignedBox> aggbucketsLEFT(k);
    std::vector<AxisAlignedBox> aggbucketsRIGHT(k);
    std::vector<int> countLEFT(k);
    std::vector<int> countRIGHT(k);

    std::sort(primitives.begin(), primitives.end(), [&](const Primitive& a, const Primitive& b) {
        return computePrimitiveCentroid(a)[axis] < computePrimitiveCentroid(b)[axis];
    });

    float final_cost = FLT_MAX;
    float final_index = -1;

    for (int i = 0; i < primitives.size(); ++i) {
        int index = (int)(((computePrimitiveCentroid(primitives[i])[axis] - aabb.lower[axis]) * k * (1 - FLT_EPSILON)) / axis_length);
        buckets[index].count++;
        buckets[index].aabb.lower = glm::min(buckets[index].aabb.lower, computePrimitiveAABB(primitives[i]).lower);
        buckets[index].aabb.upper = glm::max(buckets[index].aabb.upper, computePrimitiveAABB(primitives[i]).upper);
    }

    for (int i = 0; i < k; ++i) {

        if (i == 0) {
                aggbucketsLEFT[0] = buckets[0].aabb;
                countLEFT[0] = buckets[0].count;
        } else {
                aggbucketsLEFT[i] = AxisAlignedBox(glm::min(aggbucketsLEFT[i - 1].lower, buckets[i].aabb.lower),
                    glm::max(aggbucketsLEFT[i - 1].upper, buckets[i].aabb.upper));
                countLEFT[i] = countLEFT[i - 1] + buckets[i].count;
        }
    }

    for (int i = k - 1; i >= 0; --i) {
        if (i == k - 1) {
                aggbucketsRIGHT[k - 1] = buckets[k - 1].aabb;
                countRIGHT[k - 1] = buckets[k - 1].count;
        } else {
                aggbucketsRIGHT[i] = AxisAlignedBox(glm::min(aggbucketsRIGHT[i + 1].lower, buckets[i].aabb.lower),
                    glm::max(aggbucketsRIGHT[i + 1].upper, buckets[i].aabb.upper));
                countRIGHT[i] = countRIGHT[i + 1] + buckets[i].count;
        }
    }

    for (int i = 0; i < k - 1; ++i) {

        float cost;
        if (countLEFT[i] == 0 || countRIGHT[i + 1] == 0)
                continue;
        cost = countLEFT[i] * surfaceArea(aggbucketsLEFT[i]) + countRIGHT[i + 1] * surfaceArea(aggbucketsRIGHT[i + 1]);

        if (cost < final_cost) {
                final_cost = cost;
                final_index = countLEFT[i];
        }
    }

    if (final_cost == FLT_MAX || final_index == -1) {
        return splitPrimitivesByMedian(aabb, axis, primitives);
    }

    return final_index;
}