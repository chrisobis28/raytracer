#include "bvh.h"
#include "draw.h"
#include "extra.h"
#include "intersect.h"
#include "render.h"
#include "scene.h"
#include <bit>
#include <span>
#include <stack>

// *
// *
// MOTION BLUR BVH EXTENSION
// *
// Disabled by default, change the value of flag "ENABLE_MOTION_BLUR" in mesh.h to include motion blur
// source code in the compilation.
// *
// *
#if (ENABLE_MOTION_BLUR)

bool intersectRayWithBVHTime(RenderState& state, const BVHInterface& bvh, Ray& ray, HitInfo& hitInfo)
{
    // Relevant data in the constructed BVH
    std::span<const BVHInterface::Node> nodes = bvh.nodes();
    std::span<const BVHInterface::Primitive> primitives = bvh.primitives();

    // Return value
    bool is_hit = false;

    if (state.features.enableAccelStructure) {
        float original_t = ray.t;
        float final_t = std::numeric_limits<float>::infinity();
        std::stack<BVHInterface::Node> stack;
        stack.push(nodes[BVH::RootIndex]);

        while (!stack.empty()) {
            BVHInterface::Node curNode = stack.top();
            stack.pop();

            // Check if the ray is within physical and time bounds of the AABB
            if (ray.time >= curNode.aabb.time[0] && ray.time <= curNode.aabb.time[1] && intersectRayWithShape(curNode.aabb, ray)) {
                ray.t = original_t;
                if (curNode.isLeaf()) {
                    for (uint32_t i = 0; i < curNode.primitiveCount(); ++i) {
                        BVHInterface::Primitive p = primitives[curNode.primitiveOffset() + i];
                        Transform transform =  state.scene.meshes[p.meshID].transform;
                        p.v0 = applyTransform(p.v0, transform, ray.time);
                        p.v1 = applyTransform(p.v1, transform, ray.time);
                        p.v2 = applyTransform(p.v2, transform, ray.time);

                        // Test if ray intersects with a triangle projected in time
                        if (intersectRayWithTriangle(p.v0.position,p.v1.position, p.v2.position, ray, hitInfo)) {
                            is_hit = true;
                            if (final_t > ray.t) {
                                final_t = ray.t;
                                updateHitInfo(state, p, ray, hitInfo);
                            }
                        }
                        ray.t = original_t;
                    }
                } else {
                    stack.push(nodes[curNode.leftChild()]);
                    stack.push(nodes[curNode.rightChild()]);
                }
            }
        }

        if (is_hit) {
            ray.t = final_t;
        }
    } else {
        // Naive implementation; simply iterates over all primitives
        // Use the given time to transform the primitives according to their mesh transform function
        for (const auto& prim : primitives) {
            Transform transform =  state.scene.meshes[prim.meshID].transform;
            const Vertex& v0 = applyTransform(prim.v0, transform, ray.time);
            const Vertex& v1 = applyTransform(prim.v1, transform, ray.time);
            const Vertex& v2 = applyTransform(prim.v2, transform, ray.time);
            if (intersectRayWithTriangle(v0.position, v1.position, v2.position, ray, hitInfo)) {
                updateHitInfo(state, prim, ray, hitInfo);
                is_hit = true;
            }
        }
    }

    // Intersect with transformed spheres.
    for (const auto& sphere : state.scene.spheres) {
        Sphere transformedSphere { sphere };
        transformedSphere.center = applyTransform(transformedSphere.center, transformedSphere.transform, ray.time);
        is_hit |= intersectRayWithShape(transformedSphere, ray, hitInfo);
    }

    return is_hit;
}

AxisAlignedBox computeSpanAABBTime(std::span<const BVHInterface::Primitive> primitives, const Scene& scene, glm::vec2 time)
{
    //condition size is 0
    if (primitives.empty())
        return {};

    //initializing all absolute values with the first vertex of them all
    glm::vec3 aabbmin = { std::numeric_limits<float>::infinity(), std::numeric_limits<float>::infinity(), std::numeric_limits<float>::infinity() };
    glm::vec3 aabbmax = { -std::numeric_limits<float>::infinity(), -std::numeric_limits<float>::infinity(), -std::numeric_limits<float>::infinity() };

    //going through all the vertices and computing absolute minimum and maximum coordinates
    for (BVHInterface::Primitive p : primitives) {
        Transform transform =  scene.meshes[p.meshID].transform;

        aabbmax = glm::max(aabbmax, applyTransform(p.v0.position, transform, time[0]));
        aabbmax = glm::max(aabbmax, applyTransform(p.v0.position, transform, time[1]));
        aabbmax = glm::max(aabbmax, applyTransform(p.v1.position, transform, time[0]));
        aabbmax = glm::max(aabbmax, applyTransform(p.v1.position, transform, time[1]));
        aabbmax = glm::max(aabbmax, applyTransform(p.v2.position, transform, time[0]));
        aabbmax = glm::max(aabbmax, applyTransform(p.v2.position, transform, time[1]));

        aabbmin = glm::min(aabbmin, applyTransform(p.v0.position, transform, time[0]));
        aabbmin = glm::min(aabbmin, applyTransform(p.v0.position, transform, time[1]));
        aabbmin = glm::min(aabbmin, applyTransform(p.v1.position, transform, time[0]));
        aabbmin = glm::min(aabbmin, applyTransform(p.v1.position, transform, time[1]));
        aabbmin = glm::min(aabbmin, applyTransform(p.v2.position, transform, time[0]));
        aabbmin = glm::min(aabbmin, applyTransform(p.v2.position, transform, time[1]));
    }
    return {
        .lower = aabbmin,
        .upper = aabbmax,
        .time = time
    };
}

void BVH::buildRecursiveTime(const Scene& scene, const Features& features, std::span<Primitive> primitives, uint32_t nodeIndex, uint32_t level, glm::vec2 time)
{
    // Compute aabb spanning the whole time interval
    AxisAlignedBox aabbTime = computeSpanAABBTime(primitives, scene, time);
    // Compute aabb at single time instance at the start of the interval
    AxisAlignedBox aabb = computeSpanAABBTime(primitives, scene, glm::vec2 { time[0], time[0] });

    // Additional condition to prevent stack overflow: BVH must have a maximum of 13 levels
    if (primitives.size() <= BVH::LeafSize || level > 12) {
        m_nodes[nodeIndex] = buildLeafData(scene, features, aabbTime, primitives);

        if (level > m_numLevels) {
            m_numLevels = level+1;
        }

        m_numLeaves++;
    }
    else {
        // Split node based on the longest axis heuristic
        // The node will be split either in x/y/z dimension or in time

        uint32_t longestAxisTime = computeAABBLongestAxis(aabbTime);
        uint32_t longestAxis = computeAABBLongestAxis(aabb);
        float longestAxisLengthTime = aabbTime.upper[longestAxisTime] - aabbTime.lower[longestAxisTime];
        float longestAxisLength = aabb.upper[longestAxis] - aabb.lower[longestAxis];
        float costFactor = 2.5f;

        if (longestAxisLengthTime < longestAxisLength * costFactor) {
            // Split node in a physical dimension
            size_t index = splitPrimitivesByMedian(aabb, longestAxis, primitives);
            uint32_t indexLeft = nextNodeIdx();
            uint32_t indexRight = nextNodeIdx();

            m_nodes[nodeIndex] = buildNodeData(scene, features, aabbTime, indexLeft, indexRight);

            buildRecursiveTime(scene, features, primitives.subspan(0, index), indexLeft, level+1, time);
            buildRecursiveTime(scene, features, primitives.subspan(index, primitives.size() - index), indexRight, level+1, time);
        } else {
            // Split node in time
            float timeBoundary = (time[0] + time[1]) * 0.5f;
            uint32_t indexLeft = nextNodeIdx();
            uint32_t indexRight = nextNodeIdx();

            m_nodes[nodeIndex] = buildNodeData(scene, features, aabbTime, indexLeft, indexRight);

            buildRecursiveTime(scene, features, primitives, indexLeft, level+1, glm::vec2 { time[0], timeBoundary });
            buildRecursiveTime(scene, features, primitives, indexRight, level+1, glm::vec2 { timeBoundary, time[1] });
        }
    }
}

void BVH::debugDrawLevelTime(int level) {
    // Example showing how to draw an AABB as a (white) wireframe box.
    // Hint: use draw functions (see `draw.h`) to draw the contained boxes with different
    // colors, transparencies, etc.
    //AxisAlignedBox aabb { .lower = glm::vec3(0.0f), .upper = glm::vec3(0.0f, 1.05f, 1.05f) };
    //drawAABB(aabb, DrawMode::Wireframe, glm::vec3(0.05f, 1.0f, 0.05f), 0.1f);
    std::stack<int> toDraw;
    BVH::collectIndexes(BVH::RootIndex, 0, level, toDraw);

    int count = toDraw.size();
    while (!toDraw.empty()) {
        int nodeIndex = toDraw.top();
        toDraw.pop();
        float col = (float) toDraw.size() / (float) count * 0.5f;

        if (!m_nodes[nodeIndex].isLeaf() && m_nodes[m_nodes[nodeIndex].leftChild()].aabb.time[1] != m_nodes[nodeIndex].aabb.time[1]) {
            // This is a temporal split node
            drawAABB(m_nodes[nodeIndex].aabb, DrawMode::Wireframe, glm::vec3(0.5f + col, 0.0f, 0.0f), 1.0f );
        } else {
            // This is a spatial split node
            drawAABB(m_nodes[nodeIndex].aabb, DrawMode::Wireframe, glm::vec3(0.0f, 0.0f, 0.5f + col), 1.0f );
        }

    }
}

void BVH::collectIndexes(int nodeIndex, int currentLevel, int level, std::stack<int>& toDraw) {
    if (currentLevel == level) {
        toDraw.push(nodeIndex);
    } else {
        if (!m_nodes[nodeIndex].isLeaf()) {
            collectIndexes(m_nodes[nodeIndex].leftChild(), currentLevel + 1, level, toDraw);
            collectIndexes(m_nodes[nodeIndex].rightChild(), currentLevel + 1, level, toDraw);
        }
    }
}

#endif
