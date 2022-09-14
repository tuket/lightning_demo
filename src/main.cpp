#include <stdio.h>
#include <sokol_app.h>
#include <sokol_gfx.h>
#include <sokol_glue.h>

/* application state */
static struct {
    sg_pipeline pip;
    sg_bindings bind;
    sg_pass_action pass_action;
} state;

static void init(void) {

    const sg_desc sgDesc{
        .context = sapp_sgcontext()
    };
    sg_setup(&sgDesc);
    //__dbgui_setup(sapp_sample_count());

    /* a vertex buffer with 3 vertices */
    float vertices[] = {
        // positions            // colors
         0.0f,  0.5f, 0.5f,     1.0f, 0.0f, 0.0f, 1.0f,
         0.5f, -0.5f, 0.5f,     0.0f, 1.0f, 0.0f, 1.0f,
        -0.5f, -0.5f, 0.5f,     0.0f, 0.0f, 1.0f, 1.0f
    };
    const sg_buffer_desc bufferDesc{
        .data = SG_RANGE(vertices),
        .label = "triangle-vertices"
    };
    state.bind.vertex_buffers[0] = sg_make_buffer(&bufferDesc);

    /* create shader from code-generated sg_shader_desc */
    sg_shader shd = sg_make_shader(triangle_shader_desc(sg_query_backend()));

    /* create a pipeline object (default render states are fine for triangle) */
    const sg_pipeline_desc pipelineDesc{
        .shader = shd,
        /* if the vertex layout doesn't have gaps, don't need to provide strides and offsets */
        .layout = {
            .attrs = {
                [ATTR_vs_position].format = SG_VERTEXFORMAT_FLOAT3,
                [ATTR_vs_color0].format = SG_VERTEXFORMAT_FLOAT4
            }
        },
        .label = "triangle-pipeline"
    };
    state.pip = sg_make_pipeline(pipelineDesc);

    /* a pass action to clear framebuffer to black */
    state.pass_action = sg_pass_action{
        .colors = {{.action = SG_ACTION_CLEAR, .value = {0.0f, 0.0f, 0.0f, 1.0f } }}
    };
}

void frame(void) {
    sg_begin_default_pass(&state.pass_action, sapp_width(), sapp_height());
    sg_apply_pipeline(state.pip);
    sg_apply_bindings(&state.bind);
    sg_draw(0, 3, 1);
    //__dbgui_draw();
    sg_end_pass();
    sg_commit();
}

void cleanup(void) {
    //__dbgui_shutdown();
    sg_shutdown();
}

sapp_desc sokol_main(int argc, char* argv[]) {
    (void)argc; (void)argv;
    return sapp_desc {
        .init_cb = init,
        .frame_cb = frame,
        .cleanup_cb = cleanup,
        //.event_cb = __dbgui_event,
        .width = 640,
        .height = 480,
        .window_title = "Triangle (sokol-app)",
        .icon = {.sokol_default = true },
        .gl_force_gles2 = true,
    };
}