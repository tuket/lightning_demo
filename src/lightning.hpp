#include <glm/glm.hpp>
#include <vector>
#include <span>

struct LightningVert2d {
    glm::vec2 pos;
    float power;
};
struct LightningVert3d {
    glm::vec3 pos;
    glm::vec4 color;
};

struct LightningParams {
    float height = 0;
    float chaos = 0.3;
    float ramificationProbability = 0.45;
    glm::vec2 ramificationChaos = { 0.7, 1.5 };
    glm::vec2 ramificationLength = { 0.4, 0.8 };
    glm::vec2 ramificationPower = { 0.1, 0.3 };
    float maxThickness = 0.01;
    int maxSubdivs = 3;
    glm::vec3 cameraPos = { NAN, NAN, NAN }; // if camera position is provided, we orient the lightning geometry towards the camera (billboarding)
};

struct Lightning_linesToTriangles_inParams {
    std::span<const LightningVert2d> verts;
    std::span<const uint32_t> inds;
    LightningParams params;
};
struct Lightning_linesToTriangles_outParams {
    std::vector<LightningVert3d>& verts;
    std::vector<uint32_t>& inds;
};

void generateLighting2d(std::vector<LightningVert2d>& verts, std::vector<uint32_t>& inds,
    const LightningParams& params,
    const glm::vec2& ori, const glm::vec2& dst,
    float power = 1, int level = 0, int parent = -1);

void generateLightning_triangles(
    std::vector<LightningVert3d>& verts, std::vector<uint32_t>& inds,
    const LightningParams& params,
    const glm::vec2& ori, const glm::vec2& dst
);

void generateLightning_lines(
    std::vector<LightningVert3d>& verts, std::vector<uint32_t>& inds,
    const LightningParams& params,
    const glm::vec2& ori, const glm::vec2& dst
);

void lightning_linesToTriangles(
    Lightning_linesToTriangles_outParams& out,
    const Lightning_linesToTriangles_inParams& in);

void lightning_linesTo3d(
    std::vector<LightningVert3d>& out,
    std::span<const LightningVert2d> in);