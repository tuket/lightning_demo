@vs vs
in vec3 a_pos;
in vec4 a_color;

out vec4 v_color;

uniform vs_params {
    mat4 u_modelViewProj;
};

void main()
{
    gl_Position = u_modelViewProj * vec4(a_pos, 1);
    v_color = a_color;
}
@end

@fs fs
out vec4 o_color;

in vec4 v_color;

uniform fs_params {
    vec4 u_color;
};

void main()
{
    o_color = v_color * u_color;
}
@end

@program triangle vs fs