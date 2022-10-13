#include <stdio.h>
#define SOKOL_IMPL
#ifdef __EMSCRIPTEN__
    #define SOKOL_GLES2
#else
    #define SOKOL_GLCORE33
#endif
#include <sokol_app.h>
#include <sokol_time.h>
#include <sokol_gfx.h>
#include <sokol_glue.h>
#include <imgui.h>
#include <sokol_gfx_imgui.h>
#include <sokol_imgui.h>
#include "lightning.hpp"
#include "lightning.shader.h"
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <vector>
#include <span>
#include <stdlib.h>

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
LightningParams params;

// application state
static struct {
    sg_pipeline lines_pipeline, triangles_pipeline, niLines_pipeline; // niLines: non-indexed lines
    sg_bindings lines_bind, triangles_bind;
    sg_pass_action pass_action;
    int numLightnings = 1;
    size_t lines_numVerts = 0, lines_numInds = 0;
    size_t triangles_numVerts = 0, triangles_numInds = 0;
    std::vector<vec2> targetPoints;
    float lightningUpdatePeriod = 0.4;
    bool needToUpdateLightnings = true;
    bool needToUpdateLightningsMesh = true;
} state;

static bool watchFromRight = false;
static bool polygonal = true;

static float randFloat()
{
    return float(rand()) / RAND_MAX;
}

static vec2 randomDir2()
{
    const float a = 2.f * PI * randFloat();
    return { cos(a), sin(a) };
}

static void generateRandomLightningDests()
{
    state.targetPoints.clear();
    state.targetPoints.push_back(vec2{0.99, 0});
    for (int i = 0; i < state.numLightnings-1; i++)
        state.targetPoints.push_back(randomDir2() * glm::mix(0.5f, 1.0f, randFloat()));
}

static void init(void)
{
    stm_setup();
    generateRandomLightningDests();

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
        .data = {nullptr, numVertsEstimation * sizeof(LightningVert3d)},
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
        .data = {nullptr, 3 * numVertsEstimation * sizeof(LightningVert3d)},
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

    // create shader from code-generated sg_shader_desc
    sg_shader shd = sg_make_shader(triangle_shader_desc(sg_query_backend()));

    auto setCommonPipelineLayout = [](sg_layout_desc& layout) {
        layout.buffers[0].stride = sizeof(LightningVert3d);
        layout.attrs[ATTR_vs_a_pos] = {
            .buffer_index = 0,
            .offset = offsetof(LightningVert3d, pos),
            .format = SG_VERTEXFORMAT_FLOAT3,
        };
        layout.attrs[ATTR_vs_a_color] = {
            .buffer_index = 0,
            .offset = offsetof(LightningVert3d, color),
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

    // a pass action to clear framebuffer to black
    state.pass_action = sg_pass_action{
        .colors = {{.action = SG_ACTION_CLEAR, .value = {0.0f, 0.0f, 0.0f, 1.0f } }}
    };
}

static void updateLightnings(bool updateLines)
{
    static std::vector<LightningVert2d> lines_verts2d;
    static std::vector<LightningVert3d> lines_verts;
    static std::vector<u32> lines_inds;
    if (updateLines) {
        lines_verts2d.clear();
        lines_verts.clear();
        lines_inds.clear();
        for (const vec2& p : state.targetPoints)
            generateLighting2d(lines_verts2d, lines_inds, params, vec2(0), p);
        lightning_linesTo3d(lines_verts, lines_verts2d);

        sg_update_buffer(state.lines_bind.vertex_buffers[0], { lines_verts.data(), lines_verts.size() * sizeof(LightningVert3d) });
        state.lines_numVerts = lines_verts.size();
        sg_update_buffer(state.lines_bind.index_buffer, { lines_inds.data(), lines_inds.size() * sizeof(u32) });
        state.lines_numInds = lines_inds.size();
    }

    static std::vector<LightningVert3d> triangles_verts;
    static std::vector<u32> triangles_inds;
    triangles_verts.clear();
    triangles_inds.clear();
    Lightning_linesToTriangles_outParams outParams{ triangles_verts, triangles_inds };
    const Lightning_linesToTriangles_inParams inParams{
        .verts = lines_verts2d,
        .inds = lines_inds,
        .params = params
    };
    lightning_linesToTriangles(outParams, inParams);

    sg_update_buffer(state.triangles_bind.vertex_buffers[0], { triangles_verts.data(), triangles_verts.size() * sizeof(LightningVert3d) });
    state.triangles_numVerts = triangles_verts.size();
    sg_update_buffer(state.triangles_bind.index_buffer, { triangles_inds.data(), triangles_inds.size() * sizeof(u32) });
    state.triangles_numInds = triangles_inds.size();
}

static void drawGui()
{
    ImGui::SetNextWindowPos({ 0, 0 }, ImGuiCond_Once);
    ImGui::Begin("params", 0, ImGuiWindowFlags_AlwaysAutoResize);

    ImGui::Checkbox("watch from right", &watchFromRight);
    state.needToUpdateLightnings |= ImGui::Button("update");
    ImGui::SameLine();
    state.needToUpdateLightningsMesh |= ImGui::Button("recalc mesh");
    const int maxLightnings = 32;
    if (ImGui::SliderInt("numLightnings", &state.numLightnings, 1, maxLightnings)) {
        state.numLightnings = glm::clamp(state.numLightnings, 1, maxLightnings);
        generateRandomLightningDests();
        state.needToUpdateLightnings = true;
    }
    state.needToUpdateLightnings |= ImGui::DragFloat("update period", &state.lightningUpdatePeriod, 0.01f, 0.01f, 100);
    state.needToUpdateLightnings |= ImGui::SliderFloat("chaos factor", &params.chaos, 0, 1);
    state.needToUpdateLightnings |= ImGui::SliderFloat("ramification probability", &params.ramificationProbability, 0, 1);
    if (ImGui::SliderFloat2("ramification chaos", &params.ramificationChaos[0], 0, 3)) {
        params.ramificationChaos[0] = glm::clamp(params.ramificationChaos[0], 0.f, params.ramificationChaos[1]);
        params.ramificationChaos[1] = glm::clamp(params.ramificationChaos[1], params.ramificationChaos[1], 3.f);
        state.needToUpdateLightnings = true;
    }
    if (ImGui::SliderFloat2("ramification length", &params.ramificationLength[0], 0, 1)) {
        params.ramificationLength[0] = glm::clamp(params.ramificationLength[0], 0.f, params.ramificationLength[1]);
        params.ramificationLength[1] = glm::clamp(params.ramificationLength[1], params.ramificationLength[1], 1.f);
        state.needToUpdateLightnings = true;
    }
    if (ImGui::SliderFloat2("ramification power", &params.ramificationPower[0], 0, 1)) {
        params.ramificationPower[0] = glm::clamp(params.ramificationPower[0], 0.f, params.ramificationPower[1]);
        params.ramificationPower[1] = glm::clamp(params.ramificationPower[1], params.ramificationPower[1], 1.f);
        state.needToUpdateLightnings = true;
    }

    bool orientToCam = !isnan(params.cameraPos.x);
    if(ImGui::Checkbox("orient to cam", &orientToCam)) {
        params.cameraPos = orientToCam ? camera.pos : vec3(NAN);
        state.needToUpdateLightningsMesh = true;
    }

    ImGui::Checkbox("polygonal", &polygonal);
    ImGui::Indent();
    state.needToUpdateLightnings |= ImGui::DragFloat("maxThickness", &params.maxThickness, 0.001f);

    ImGui::End();
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
    if (polygonal) {
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

    if (polygonal)
        sg_draw(0, state.triangles_numInds, 1);
    else
        sg_draw(0, state.lines_numInds, 1);

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