#include "interpolate.h"
#include <glm/geometric.hpp>

// TODO Standard feature
// Given three triangle vertices and a point on the triangle, compute the corresponding barycentric coordinates of the point.
// and return a vec3 with the barycentric coordinates (alpha, beta, gamma).
// - v0;     Triangle vertex 0
// - v1;     Triangle vertex 1
// - v2;     Triangle vertex 2
// - p;      Point on triangle
// - return; Corresponding barycentric coordinates for point p.
// This method is unit-tested, so do not change the function signature.
glm::vec3 computeBarycentricCoord(const glm::vec3& v0, const glm::vec3& v1, const glm::vec3& v2, const glm::vec3& p)
{
    glm::vec3 u0 = v0 - v2;
    glm::vec3 u1 = v1 - v2;
    glm::vec3 u2 = p - v2;

    // The equation for finding barycentric coordinates is
    // a * u0 + b * u1 = u2
    // if x = [a, b] and M = [u0, u1] then M*x = u2
    // multiplying both sides by Mt (M transposed), we have Mt*M*x = Mt*u2, a 2x2 linear system, which can be solved using Cramer's rule
    //
    // *** The derivation of this linear system is based on a book by Christer Ericson's: Real-Time Collision Detection ***
    float m00 = glm::dot(u0, u0);
    float m01 = glm::dot(u0, u1);
    float m02 = glm::dot(u0, u2);
    float m10 = m01;
    float m11 = glm::dot(u1, u1);
    float m12 = glm::dot(u1, u2);

    float det = m00 * m11 - m10 * m01;
    float det_a = m02 * m11 - m12 * m01;
    float det_b = m00 * m12 - m01 * m02;

    float a = det_a / det;
    float b = det_b / det;

    return glm::vec3 { a, b, 1.f - a - b };
}

// TODO Standard feature
// Linearly interpolate three normals using barycentric coordinates.
// - n0;     Triangle normal 0
// - n1;     Triangle normal 1
// - n2;     Triangle normal 2
// - bc;     Barycentric coordinate
// - return; The smoothly interpolated normal.
// This method is unit-tested, so do not change the function signature.
glm::vec3 interpolateNormal(const glm::vec3& n0, const glm::vec3& n1, const glm::vec3& n2, const glm::vec3 bc)
{
    return bc.x * n0 + bc.y * n1 + bc.z * n2;
}

// TODO Standard feature
// Linearly interpolate three texture coordinates using barycentric coordinates.
// - n0;     Triangle texture coordinate 0
// - n1;     Triangle texture coordinate 1
// - n2;     Triangle texture coordinate 2
// - bc;     Barycentric coordinate
// - return; The smoothly interpolated texturre coordinate.
// This method is unit-tested, so do not change the function signature.
glm::vec2 interpolateTexCoord(const glm::vec2& t0, const glm::vec2& t1, const glm::vec2& t2, const glm::vec3 bc)
{
    return bc.x * t0 + bc.y * t1 + bc.z * t2;
}
