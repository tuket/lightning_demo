#include <stdio.h>
#define SOKOL_IMPL
#define SOKOL_GLCORE33
#include <sokol_app.h>
#include <sokol_time.h>
#include <sokol_gfx.h>
#include <sokol_glue.h>
#include <imgui.h>
#include <sokol_gfx_imgui.h>
#include <sokol_imgui.h>
#include "lightning.shader.h"
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <vector>
#include <span>
#include <unordered_map>

#define DEBUG_LINES 0

using glm::vec2;
using glm::vec3;
using glm::vec4;
constexpr static float PI = 3.141592653589;
typedef uint32_t u32;
typedef uint64_t u64;

static const vec4 RED = { 1, 0, 0, 1 };
static const vec4 GREEN = { 0, 1, 0, 1 };
static const vec4 BLUE = { 0, 0, 1, 1 };
static const vec4 WHITE = { 1, 1, 1, 1 };
static const vec4 k_transparent = { 1, 1, 1, 0 };
static const vec4 k_opaque = { 1,1,1, 1 };

struct Vert3d {
    vec3 pos;
    vec4 color;
};

struct LightningVert2d {
    vec2 pos;
    float power;
};

struct Camera {
    vec3 pos;
    vec3 target;
    float fovY;
    float zNear, zFar;
};
static const Camera camera = {
    .pos = {1, 1, 1},
    .target = {0, 0, 0},
    .fovY = glm::radians(60.f),
    .zNear = 0.1, .zFar = 100,
};

/* application state */
static struct {
    sg_pipeline lines_pipeline, triangles_pipeline, niLines_pipeline; // niLines: non-indexed lines
    sg_bindings lines_bind, triangles_bind;
#if DEBUG_LINES
    sg_bindings debugLines_bind;
#endif
    sg_pass_action pass_action;
    size_t lines_numVerts = 0, lines_numInds = 0;
    size_t triangles_numVerts = 0, triangles_numInds = 0;
    std::vector<vec2> targetPoints;
    float lightningUpdatePeriod = 1000;
    bool needToUpdateLightnings = true;
    bool needToUpdateLightningsMesh = false;
} state;

struct LightningParams {
    float chaos = 0.2;
    //float ramificationProbability = 0.25;
    float ramificationProbability = 0.0;
    float ramificationChaos[2] = { 0.1, 0.7 };
    float ramificationLength[2] = { 0.2, 0.8 };
    float ramificationPower[2] = { 0.1, 0.3 };
    float maxThickness = 0.01;
    vec3 cameraPos = {NAN, NAN, NAN}; // if camera position is provided, we orient the lightning geometry towards the camera (billboarding)
    int maxSubdivs = 3;
    bool polygonal = false;
} params;

struct DebugLine {
    vec2 ori, dir;
    vec4 color;
};

static bool watchFromRight = false;
#if DEBUG_LINES
static std::vector<std::vector<DebugLine>> debugLines;
static size_t selectedDebugLine = -1;
#endif

static float randFloat()
{
    return float(rand()) / RAND_MAX;
}

static vec2 randomDir2()
{
    const float a = 2.f * PI * randFloat();
    return { cos(a), sin(a) };
}

static vec3 randomDir3()
{
    const float y = -1 + 2.f * randFloat();
    const float phi = 2.f * PI * randFloat();
    const float x = cos(phi);
    const float z = sin(phi);
    const float xy = sqrtf(1 - y*y);
    return { x * xy, y, z * xy };
}

static vec3 randomOrthoDir(vec3 axis)
{
    vec3 X;
    if (dot(axis, vec3(1, 0, 0)) > 0.75)
        X = normalize(cross(axis, vec3(1, 0, 0)));
    else
        X = normalize(cross(axis, vec3(0, 1, 0)));

    const vec3 Y = cross(axis, X);
    const vec2 R = randomDir2();
    return R.x * X + R.y * Y;
}

/*static void createLighting3d(std::vector<LightningVert3d>& verts, std::vector<u32>& inds,
    const LightningParams& params,
    const vec3& ori, const vec3& dst, float power = 1, int level = 0, int parent = -1)
{
    if (level >= params.maxSubdivs) {
        inds.push_back(parent);
        inds.push_back(verts.size());
        verts.push_back({dst, power});
        return;
    }
    if (level == 0) {
        verts.push_back({ ori, power });
        parent = 0;
    }

    const vec3 od = dst - ori;
    const float odLen = length(od);
    const vec3 dispDir = randomOrthoDir(od / odLen);

    const float midSplitPercent = 0.4f + 0.2f * randFloat();
    vec3 midP = glm::mix(ori, dst, midSplitPercent);
    midP += dispDir * odLen * params.chaos;

    createLighting3d(verts, inds, params, ori, midP, power, level + 1, parent);
    const size_t midP_ind = verts.size() - 1;

    if (randFloat() < params.ramificationProbability) {
        vec3 ramiDir = 2 * randFloat() * params.ramificationLength * (dst - midP);
        const float ramiDirLen = length(ramiDir);
        ramiDir += params.chaos * ramiDirLen * randomOrthoDir(ramiDir / ramiDirLen);
        vec3 ramiDst = midP + ramiDir;
        createLighting3d(verts, inds, params, midP, ramiDst, power * 0.3f, level + 1, midP_ind);
    }

    createLighting3d(verts, inds, params, midP, dst, power, level + 1, midP_ind);
}*/

static void createLighting2d(std::vector<LightningVert2d>& verts, std::vector<u32>& inds,
    const LightningParams& params,
    const vec2& ori, const vec2& dst, float power = 1, int level = 0, int parent = -1)
{
    if (level >= params.maxSubdivs) {
        inds.push_back(parent);
        inds.push_back(verts.size());
        verts.push_back({ dst, power });
        return;
    }
    if (level == 0) {
        verts.push_back({ ori, power });
        parent = 0;
    }

    vec2 od = dst - ori;
    const float odLen = length(od);
    od /= odLen;
    const vec2 dispDir = {-od.y, od.x};

    const float midSplitPercent = 0.4f + 0.2f * randFloat();
    vec2 midP = glm::mix(ori, dst, midSplitPercent);
    midP += (randFloat() - 0.5f) * dispDir * odLen * params.chaos;

    createLighting2d(verts, inds, params, ori, midP, power, level + 1, parent);
    const size_t midP_ind = verts.size() - 1;

    if (randFloat() < params.ramificationProbability) {
        vec2 ramiDir = glm::mix(params.ramificationLength[0], params.ramificationLength[1], randFloat()) * power * (dst - midP);
        const float ramiDirLen = length(ramiDir);
        ramiDir /= ramiDirLen;
        const float randChaos = 2 * randFloat() - 1;
        const float ramiChaos = glm::sign(randChaos) * glm::mix(params.ramificationChaos[0], params.ramificationChaos[1], abs(randChaos));
        ramiDir += ramiDir + ramiChaos * vec2(-ramiDir.y, ramiDir.x);
        ramiDir = normalize(ramiDir);
        ramiDir *= ramiDirLen;
        vec2 ramiDst = midP + ramiDir;
        const float ramiPower = glm::mix(params.ramificationPower[0], params.ramificationPower[1], randFloat()) * power;
        createLighting2d(verts, inds, params, midP, ramiDst, ramiPower, level + 1, midP_ind);
        power -= ramiPower;
    }

    createLighting2d(verts, inds, params, midP, dst, power, level + 1, midP_ind);
}

static void init(void)
{
    stm_setup();
    for (int i = 0; i < 1; i++)
        state.targetPoints.push_back(randomDir3());

    const sg_desc sgDesc{
        .context = sapp_sgcontext()
    };
    sg_setup(&sgDesc);
    const simgui_desc_t simguiDesc = {

    };
    simgui_setup(&simguiDesc);

    // lightning lines buffer
    const u32 numVertsEstimation = 16 << 10;
    const sg_buffer_desc lines_vertBufferDesc{
        .usage = SG_USAGE_STREAM,
        .data = {nullptr, numVertsEstimation * sizeof(Vert3d)},
        .label = "lines_vertices",
    };
    state.lines_bind.vertex_buffers[0] = sg_make_buffer(&lines_vertBufferDesc);

    const sg_buffer_desc lines_indBufferDesc{
        .type = SG_BUFFERTYPE_INDEXBUFFER,
        .usage = SG_USAGE_STREAM,
        .data = {nullptr, numVertsEstimation * sizeof(u32)},
        .label = "lines_indices",
    };
    state.lines_bind.index_buffer = sg_make_buffer(&lines_indBufferDesc);

    // lightning triangles buffer
    const sg_buffer_desc triangles_vertBufferDesc{
        .usage = SG_USAGE_STREAM,
        .data = {nullptr, 3 * numVertsEstimation * sizeof(Vert3d)},
        .label = "triangles_vertices",
    };
    state.triangles_bind.vertex_buffers[0] = sg_make_buffer(&triangles_vertBufferDesc);

    const sg_buffer_desc triangles_indBufferDesc{
        .type = SG_BUFFERTYPE_INDEXBUFFER,
        .usage = SG_USAGE_STREAM,
        .data = {nullptr, 3 * numVertsEstimation * sizeof(u32)},
        .label = "triangles_indices",
    };
    state.triangles_bind.index_buffer = sg_make_buffer(&triangles_indBufferDesc);

    // debug lines buffer
#if DEBUG_LINES
    const sg_buffer_desc debugLines_vertBufferDesc{
        .usage = SG_USAGE_STREAM,
        .data = {nullptr, 6 * numVertsEstimation * sizeof(Vert3d)},
        .label = "debugLines_vertices",
    };
    state.debugLines_bind.vertex_buffers[0] = sg_make_buffer(&lines_vertBufferDesc);
#endif

    /* create shader from code-generated sg_shader_desc */
    sg_shader shd = sg_make_shader(triangle_shader_desc(sg_query_backend()));

    auto setCommonPipelineLayout = [](sg_layout_desc& layout) {
        layout.buffers[0].stride = sizeof(Vert3d);
        layout.attrs[ATTR_vs_a_pos] = {
            .buffer_index = 0,
            .offset = offsetof(Vert3d, pos),
            .format = SG_VERTEXFORMAT_FLOAT3,
        };
        layout.attrs[ATTR_vs_a_color] = {
            .buffer_index = 0,
            .offset = offsetof(Vert3d, color),
            .format = SG_VERTEXFORMAT_FLOAT4,
        };
    };

    // lines pipelines
    sg_pipeline_desc lines_pipelineDesc {
        .shader = shd,
        .primitive_type = SG_PRIMITIVETYPE_LINES,
        .index_type = SG_INDEXTYPE_UINT32,
        /* if the vertex layout doesn't have gaps, don't need to provide strides and offsets */
        .label = "lines_pipeline",
    };
    setCommonPipelineLayout(lines_pipelineDesc.layout);
    state.lines_pipeline = sg_make_pipeline(lines_pipelineDesc);

    // triangles pipeline
    sg_pipeline_desc triangles_pipelineDesc {
        .shader = shd,
        .primitive_type = SG_PRIMITIVETYPE_TRIANGLES,
        .index_type = SG_INDEXTYPE_UINT32,
        .label = "triangles_pipeline",
    };
    triangles_pipelineDesc.colors[0].blend = {
        .enabled = true,
        .src_factor_rgb = SG_BLENDFACTOR_SRC_ALPHA,
        .dst_factor_rgb = SG_BLENDFACTOR_ONE_MINUS_SRC_ALPHA,
        .op_rgb = SG_BLENDOP_ADD,
    };
    setCommonPipelineLayout(triangles_pipelineDesc.layout);
    state.triangles_pipeline = sg_make_pipeline(triangles_pipelineDesc);

    // non-indexed lines pipeline
    sg_pipeline_desc niLines_pipelineDesc{
        .shader = shd,
        .primitive_type = SG_PRIMITIVETYPE_LINES,
        .index_type = SG_INDEXTYPE_NONE,
        .label = "niLines_pipeline",
    };
    setCommonPipelineLayout(niLines_pipelineDesc.layout);
    state.niLines_pipeline = sg_make_pipeline(niLines_pipelineDesc);

    /* a pass action to clear framebuffer to black */
    state.pass_action = sg_pass_action{
        .colors = {{.action = SG_ACTION_CLEAR, .value = {0.0f, 0.0f, 0.0f, 1.0f } }}
    };
}

static vec2 to2d(vec3 p)
{
    return { p.x, -p.z };
}

static vec3 to3d(vec2 p)
{
    return {p.x, 0, -p.y};
}


static float rayVsLine(vec2 rayOri, vec2 rayDir, vec2 lineOri, vec2 lineDir)
{
    const vec2 D = lineOri - rayOri;
    const float numer = lineDir.y * D.x - lineDir.x * D.y;
    const float denom = rayDir.x * lineDir.y - rayDir.y * lineDir.x;
    return numer / denom;
}

static void lightningTo3d(std::vector<Vert3d>& out, std::span<const LightningVert2d> in)
{
    for (auto v : in) {
        out.push_back({ to3d(v.pos), WHITE });
    }
}

struct Lightning_linesToTriangles_inParams {
    std::span<const LightningVert2d> verts;
    std::span<const u32> inds;
    const vec3 camPos;
    float maxThickness;
};
struct Lightning_linesToTriangles_outParams {
    std::vector<Vert3d>& verts;
    std::vector<u32>& inds;
};

static void lightning_linesToTriangles_rec(
    Lightning_linesToTriangles_outParams& out,
    std::vector<u32>& toReorient,
    const Lightning_linesToTriangles_inParams& in,
    const std::unordered_multimap<u32, u32>& m,
    u32 parent, float parentR, u32 node, u32 parentI0, u32 parentI1, u32 parentI2)
{
    auto doQuad = [&out](u32 i0, u32 i1, u32 i2, u32 i3)
    {
        out.inds.push_back(i0);
        out.inds.push_back(i1);
        out.inds.push_back(i2);
        out.inds.push_back(i0);
        out.inds.push_back(i2);
        out.inds.push_back(i3);
    };

    const size_t numChildren = m.count(node);
    if (numChildren == 0) {
        const u32 i0 = out.verts.size();
        out.verts.push_back({ to3d(in.verts[node].pos), k_transparent});
        out.inds.push_back(i0);
        out.inds.push_back(parentI1);
        out.inds.push_back(parentI0);
        out.inds.push_back(i0);
        out.inds.push_back(parentI2);
        out.inds.push_back(parentI1);
    }
    else {
        auto childrenBegin = m.lower_bound(node);
        const float& power = in.verts[node].power;
        const float r = in.maxThickness * sqrt(power);
        const vec2& a = in.verts[parent].pos;
        const vec2& b = in.verts[node].pos;
        vec2 ab = b - a;
        const float abLen = length(ab);
        ab /= abLen;
        //const vec3 toCam = normalize(in.camPos - to3d(b));
        //const vec3 camRight = { toCam.z, 0, -toCam.x };
        //const vec3 ba_TtoCam = normalize(toCam * dot(toCam, ab) - ab);
        if (numChildren == 1) {
            const u32 child = childrenBegin->second;
            const vec2& c = in.verts[child].pos;
            const vec2 bc = normalize(c - b);
            const vec2 abRight = r * vec2(-ab.y, ab.x);
            const vec2 bcRight = r * vec2(-bc.y, bc.x);
            vec2 leftP, rightP;
            if (dot(ab, bc) < 0.99) {
                const float leftD = rayVsLine(
                    a - abRight, ab,
                    b - bcRight, bc
                );
                leftP = a - abRight + ab * leftD;
                const float rightD = rayVsLine(
                    a + abRight, ab,
                    b + bcRight, bc
                );
                rightP = a + abRight + ab * rightD;
            }
            else {
                leftP = b - abRight;
                rightP = b + abRight;
            }

            const u32 i0 = out.verts.size();
            out.verts.push_back({ to3d(leftP), k_transparent });
            out.verts.push_back({ to3d(b), k_opaque });
            out.verts.push_back({ to3d(rightP), k_transparent });

            if (!isnan(params.cameraPos.x)) {
                toReorient.push_back(i0 + 1);
                toReorient.push_back(i0);
                toReorient.push_back(i0 + 1);
                toReorient.push_back(i0 + 2);
            }

            doQuad(i0, i0 + 1, parentI1, parentI0);
            doQuad(i0 + 1, i0 + 2, parentI2, parentI1);
            lightning_linesToTriangles_rec(out, toReorient, in, m, node, r, child, i0, i0 + 1, i0 + 2);
        }
        else {
            assert(m.count(node) == 2);
            
            u32 childL = childrenBegin->second;
            auto childIt = childrenBegin;
            childIt++;
            u32 childR = childIt->second;

            const vec2 abRight = {-ab.y, b.x};
            if (dot(abRight, normalize(in.verts[childL].pos - b)) > dot(abRight, normalize(in.verts[childR].pos - b)))
                std::swap(childL, childR);
            const vec2& cl = in.verts[childL].pos;
            const vec2& cr = in.verts[childR].pos;
            const float rl = in.maxThickness * sqrt(in.verts[childL].power);
            const float rr = in.maxThickness * sqrt(in.verts[childR].power);

            vec2 bcl = cl - b;
            const float bclLen = length(bcl);
            bcl /= bclLen;
            vec2 bcr = cr - b;
            const float bcrLen = length(bcr);
            bcr /= bcrLen;
            
            const vec2 bclRight = { -bcl.y, bcl.x };
            const vec2 bcrRight = { -bcr.y, bcr.x };

            std::vector<DebugLine> DL;

            vec2 leftP;
            {
                const vec2 parentL = to2d(out.verts[parentI0].pos);
                float d = rayVsLine(
                    parentL, ab,
                    b - bclRight * rl, bcl
                );
                if (isinf(d) || isnan(d))
                    leftP = parentL;
                else {
                    d = glm::clamp(d, 0.f, 1.2f*abLen);
                    leftP = parentL + d * ab;
                }
                DL.push_back(DebugLine{ parentL, ab, RED });
                DL.push_back(DebugLine{ b - bclRight * rl, bcl, RED });
            }

            vec2 rightP;
            {
                const vec2 parentR = to2d(out.verts[parentI2].pos);
                float d = rayVsLine(
                    parentR, ab,
                    b + bcrRight * rr, bcr
                );
                if (isinf(d) || isnan(d))
                    rightP = parentR;
                else {
                    d = glm::clamp(d, 0.f, 1.2f*abLen);
                    rightP = parentR + d * ab;
                }
            }

            vec2 midP;
            {
                float d = rayVsLine(
                    b + bclRight * rl, bcl,
                    b - bcrRight * rr, bcr
                );
                if (isinf(d) || isnan(d))
                    midP = b;
                else {
                    d = glm::clamp(d, 0.f, glm::min(bclLen, bcrLen));
                    midP = b + bclRight * rl + d * bcl;
                }
                //midP = b; // hack
            }

            u32 i0 = out.verts.size();
            out.verts.push_back({ to3d(leftP), k_transparent});
            out.verts.push_back({ to3d(b), k_opaque});
            out.verts.push_back({ to3d(rightP), k_transparent});
            out.verts.push_back({ to3d(midP), k_transparent});

            if (!isnan(params.cameraPos.x)) {
                toReorient.push_back(i0 + 1);
                toReorient.push_back(i0);
                toReorient.push_back(i0 + 1);
                toReorient.push_back(i0 + 2);
                toReorient.push_back(i0 + 1);
                toReorient.push_back(i0 + 3);
            }

            doQuad(i0, i0+1, parentI1, parentI0);
            doQuad(i0+1, i0+2, parentI2, parentI1);

#if DEBUG_LINES
            debugLines.push_back(std::move(DL));
#endif

            lightning_linesToTriangles_rec(out, toReorient, in, m, node, rl, childL, i0, i0 + 1, i0 + 3);
            lightning_linesToTriangles_rec(out, toReorient, in, m, node, rr, childR, i0 + 3, i0 + 1, i0 + 2);
        }
    }
}

static void lightning_linesToTriangles(
    Lightning_linesToTriangles_outParams& out,
    const Lightning_linesToTriangles_inParams& in)
{
#if DEBUG_LINES
    debugLines.resize(0);
    selectedDebugLine = -1;
#endif

    std::unordered_multimap<u32, u32> m;
    for (size_t i = 0; i < in.inds.size(); i += 2) {
        m.insert({ in.inds[i], in.inds[i + 1] });
    }
    
    const u32 root = 0;
    assert(m.count(root) == 1);

    const u32 rootChild = m.lower_bound(root)->second;
    const vec2& a = in.verts[root].pos;
    const vec2& b = in.verts[rootChild].pos;
    const vec2 ab = normalize(b - a);
    //const vec3 toCam = normalize(in.camPos - a);
    const vec2 right = normalize(vec2(-ab.y, ab.x));
    const float& power = in.verts[root].power;
    const float r = in.maxThickness * sqrt(power);
    const vec2 al = a - right * r;
    const vec2 ar = a + right * r;
    const u32 i0 = out.verts.size();
    out.verts.push_back({ to3d(al), k_transparent });
    out.verts.push_back({ to3d(a), k_opaque });
    out.verts.push_back({ to3d(ar), k_transparent });
    std::vector<u32> toReorient;
    if (!isnan(params.cameraPos.x)) {
        toReorient.reserve(in.inds.size() * 2);
        toReorient.push_back(1);
        toReorient.push_back(0);
        toReorient.push_back(1);
        toReorient.push_back(2);
    }
    lightning_linesToTriangles_rec(out, toReorient, in, m, root, r, rootChild, i0, i0+1, i0+2);

    for (size_t i = 0; i < toReorient.size(); i += 2) {
        const vec3& midP = out.verts[toReorient[i]].pos;
        vec3& p = out.verts[toReorient[i + 1]].pos;

        const vec3 toCam = normalize(params.cameraPos - midP);
        vec3 toP = p - midP;
        const float r = length(toP);
        toP = toP - toCam * dot(toCam, toP);
        toP = r * normalize(toP);
        p = midP + toP;
    }
}

static void updateLightnings(bool updateLines)
{
    static std::vector<LightningVert2d> verts2d;
    static std::vector<u32> inds;
    if (updateLines) {
        verts2d.clear();
        inds.clear();
        for (const vec2& p : state.targetPoints)
            createLighting2d(verts2d, inds, params, { 0,0 }, p);

        std::vector<Vert3d> verts3d;
        lightningTo3d(verts3d, verts2d);

        sg_update_buffer(state.lines_bind.vertex_buffers[0], { verts3d.data(), verts3d.size() * sizeof(Vert3d) });
        state.lines_numVerts = verts3d.size();
        sg_update_buffer(state.lines_bind.index_buffer, { inds.data(), inds.size() * sizeof(u32) });
        state.lines_numInds = inds.size();
    }

    static std::vector<Vert3d> triangles_verts;
    static std::vector<u32> triangles_inds;
    triangles_verts.clear();
    triangles_inds.clear();
    {
        Lightning_linesToTriangles_outParams outParams{ triangles_verts, triangles_inds };
        const Lightning_linesToTriangles_inParams inParams{
            .verts = verts2d,
            .inds = inds,
            .camPos = camera.pos,
            .maxThickness = params.maxThickness,
        };

        lightning_linesToTriangles(outParams, inParams);
    }
    sg_update_buffer(state.triangles_bind.vertex_buffers[0], { triangles_verts.data(), triangles_verts.size() * sizeof(Vert3d) });
    state.triangles_numVerts = triangles_verts.size();
    sg_update_buffer(state.triangles_bind.index_buffer, { triangles_inds.data(), triangles_inds.size() * sizeof(u32) });
    state.triangles_numInds = triangles_inds.size();

#if DEBUG_LINES
    std::vector<Vert3d> debugLinesData;
    debugLinesData.reserve(1 << 20);
    for (const auto& v : debugLines) {
        for (const DebugLine& dl : v) {
            const vec2 dst = dl.ori + dl.dir * 1000.f;
            debugLinesData.push_back({to3d(dl.ori), dl.color});
            debugLinesData.push_back({to3d(dst), dl.color});
        }
    }
    if(debugLinesData.size())
        sg_update_buffer(state.debugLines_bind.vertex_buffers[0], { debugLinesData.data(), debugLinesData.size() * sizeof(Vert3d) });
#endif
}

static void drawGui()
{
    //ImGui::ShowDemoWindow();

    ImGui::Begin("window");

    ImGui::Checkbox("watch from right", &watchFromRight);
    state.needToUpdateLightnings |= ImGui::Button("update");
    ImGui::SameLine();
    state.needToUpdateLightningsMesh |= ImGui::Button("recalc mesh");
    state.needToUpdateLightnings |= ImGui::DragFloat("update period", &state.lightningUpdatePeriod, 0.01f, 0.01f, 100);
    state.needToUpdateLightnings |= ImGui::SliderFloat("chaos factor", &params.chaos, 0, 1);
    state.needToUpdateLightnings |= ImGui::SliderFloat("ramification probability", &params.ramificationProbability, 0, 1);
    if (ImGui::SliderFloat2("ramification chaos", params.ramificationChaos, 0, 3)) {
        params.ramificationChaos[0] = glm::clamp(params.ramificationChaos[0], 0.f, params.ramificationChaos[1]);
        params.ramificationChaos[1] = glm::clamp(params.ramificationChaos[1], params.ramificationChaos[1], 3.f);
        state.needToUpdateLightnings = true;
    }
    if (ImGui::SliderFloat2("ramification length", params.ramificationLength, 0, 1)) {
        params.ramificationLength[0] = glm::clamp(params.ramificationLength[0], 0.f, params.ramificationLength[1]);
        params.ramificationLength[1] = glm::clamp(params.ramificationLength[1], params.ramificationLength[1], 1.f);
        state.needToUpdateLightnings = true;
    }
    if (ImGui::SliderFloat2("ramification power", params.ramificationPower, 0, 1)) {
        params.ramificationPower[0] = glm::clamp(params.ramificationPower[0], 0.f, params.ramificationPower[1]);
        params.ramificationPower[1] = glm::clamp(params.ramificationPower[1], params.ramificationPower[1], 1.f);
        state.needToUpdateLightnings = true;
    }

    bool orientToCam = !isnan(params.cameraPos.x);
    if(ImGui::Checkbox("orient to cam", &orientToCam)) {
        params.cameraPos = orientToCam ? camera.pos : vec3(NAN);
        state.needToUpdateLightningsMesh = true;
    }

    ImGui::Checkbox("polygonal", &params.polygonal);
    ImGui::Indent();
    state.needToUpdateLightnings |= ImGui::DragFloat("maxThickness", &params.maxThickness, 0.001f);

    ImGui::End();

#if DEBUG_LINES
    ImGui::Begin("debug lines");

    if (ImGui::BeginTable("t", 2, ImGuiTableFlags_BordersInnerV))
    {
        ImGui::TableNextRow();
        ImGui::TableSetColumnIndex(0);
        for (size_t i = 0; i < debugLines.size(); i++) {
            bool selected = selectedDebugLine == i;
            char label[8];
            snprintf(label, std::size(label), "%d", int(i));
            const bool changed = ImGui::Selectable(label, &selected, ImGuiSelectableFlags_None);
            if (changed) {
                if (selected)
                    selectedDebugLine = i;
                else
                    selectedDebugLine = -1;
            }
        }

        ImGui::TableSetColumnIndex(1);
        ImGui::Text("properties");
        ImGui::EndTable();
    }

    ImGui::End();
#endif
}

static void frame() {
    static u64 ticks = stm_now();
    const u64 elapsedTicks = stm_laptime(&ticks);
    const float t = stm_sec(ticks);
    const float dt = stm_sec(elapsedTicks);

    static float timeSinceLightningsUpdated = 999999;
    if (timeSinceLightningsUpdated > state.lightningUpdatePeriod) {
        state.needToUpdateLightnings = true;
        timeSinceLightningsUpdated = 0;
    }
    timeSinceLightningsUpdated += dt;

    const int screenW = sapp_width();
    const int screenH = sapp_height();
    const simgui_frame_desc_t simguiFrameDesc = {
        .width = screenW, .height = screenH,
        .delta_time = 1.f / 60.f,
        .dpi_scale = 1.f,
    };
    simgui_new_frame(&simguiFrameDesc);

    drawGui();

    if (state.needToUpdateLightnings) {
        updateLightnings(true);
    }
    else if (state.needToUpdateLightningsMesh) {
        updateLightnings(false);
        state.needToUpdateLightningsMesh = false;
    }
    state.needToUpdateLightnings = false;

    sg_begin_default_pass(&state.pass_action, sapp_width(), sapp_height());
    if (params.polygonal) {
        sg_apply_pipeline(state.triangles_pipeline);
        sg_apply_bindings(&state.triangles_bind);
    }
    else {
        sg_apply_pipeline(state.lines_pipeline);
        sg_apply_bindings(&state.lines_bind);
    }

    auto applyCommonUniforms = [screenW, screenH]() {
        const float aspectRatio = float(screenW) / screenH;
        vs_params_t vsParams;
        glm::mat4 viewMtx, projMtx;
        if (watchFromRight) {
            projMtx = glm::ortho(-aspectRatio, +aspectRatio, -1.f, +1.f);
            vec3 lookDir = normalize(camera.target - camera.pos);
            lookDir = { lookDir.z, 0, -lookDir.x };
            viewMtx = glm::lookAt(camera.target + lookDir, camera.target, vec3(0, 1, 0));
        }
        else {
            projMtx = glm::perspective(camera.fovY, float(screenW) / screenH, camera.zNear, camera.zFar);
            viewMtx = glm::lookAt(camera.pos, camera.target, vec3(0, 1, 0));
        }
        const auto viewProjMtx = projMtx * viewMtx;
        memcpy(vsParams.u_modelViewProj, &viewProjMtx[0][0], sizeof(glm::mat4));
        sg_apply_uniforms(SG_SHADERSTAGE_VS, SLOT_vs_params, SG_RANGE_REF(vsParams));
    };
    applyCommonUniforms();

    fs_params_t fsParams;
    const vec4 lightningColor = { 1, 1, 0, 1 };
    memcpy(fsParams.u_color, &lightningColor[0], sizeof(vec4));
    sg_apply_uniforms(SG_SHADERSTAGE_FS, SLOT_fs_params, SG_RANGE_REF(fsParams));

    if (params.polygonal)
        sg_draw(0, state.triangles_numInds, 1);
    else
        sg_draw(0, state.lines_numInds, 1);

#if DEBUG_LINES
    if (debugLines.size()) {
        sg_apply_pipeline(state.niLines_pipeline);
        sg_apply_bindings(state.debugLines_bind);
        applyCommonUniforms();
        sg_draw(0, debugLines.size(), 1);
    }
#endif

    simgui_render();

    sg_end_pass();
    sg_commit();
}

static void cleanup()
{
    simgui_shutdown();
    sg_shutdown();
}

static void onEvent(const sapp_event* ev)
{
    simgui_handle_event(ev);
}

sapp_desc sokol_main(int argc, char* argv[]) {
    (void)argc; (void)argv;
    return sapp_desc {
        .init_cb = init,
        .frame_cb = frame,
        .cleanup_cb = cleanup,
        .event_cb = onEvent,
        .width = 1280,
        .height = 960,
        .window_title = "Lightning demo",
        .icon = {.sokol_default = true },
        .gl_force_gles2 = true,
    };
}