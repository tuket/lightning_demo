#include "lightning.hpp"
#include <array>
#include <span>

typedef uint32_t u32;
using glm::vec2, glm::vec3, glm::vec4;

// helper struct to store the child indices in a (non-complete)binary tree
struct Childs : std::array<u32, 2> {
    Childs() {
        operator[](0) = operator[](1) = -1;
    }

    size_t size()const {
        return
            this->operator[](0) == u32(-1) ? 0 :
            this->operator[](1) == u32(-1) ? 1 : 2;
    }
    void push_back(u32 x) {
        const size_t s = size();
        assert(s < 2);
        this->operator[](s) = x;
    }
};

constexpr static float PI = 3.141592653589;

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
    const float xy = sqrtf(1 - y * y);
    return { x * xy, y, z * xy };
}

static glm::vec2 to2d(glm::vec3 p)
{
    return glm::vec2(p.x, -p.z);
}

static vec3 to3d(vec2 p)
{
    return { p.x, 0, -p.y };
}

static float rayVsLine(glm::vec2 rayOri, glm::vec2 rayDir, glm::vec2 lineOri, glm::vec2 lineDir)
{
    const glm::vec2 D = lineOri - rayOri;
    const float numer = lineDir.y * D.x - lineDir.x * D.y;
    const float denom = rayDir.x * lineDir.y - rayDir.y * lineDir.x;
    return numer / denom;
}

void generateLighting2d(std::vector<LightningVert2d>& verts, std::vector<u32>& inds,
    const LightningParams& params,
    const vec2& ori, const vec2& dst, float power, int level, int parent)
{
    if (level >= params.maxSubdivs) {
        inds.push_back(parent);
        inds.push_back(verts.size());
        verts.push_back({ dst, power });
        return;
    }
    if (level == 0) {
        parent = verts.size();
        verts.push_back({ ori, power });
    }

    vec2 od = dst - ori;
    const float odLen = length(od);
    od /= odLen;
    const vec2 dispDir = { -od.y, od.x };

    const float midSplitPercent = 0.4f + 0.2f * randFloat();
    vec2 midP = glm::mix(ori, dst, midSplitPercent);
    midP += (randFloat() - 0.5f) * dispDir * odLen * params.chaos;

    generateLighting2d(verts, inds, params, ori, midP, power, level + 1, parent);
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
        generateLighting2d(verts, inds, params, midP, ramiDst, ramiPower, level + 1, midP_ind);
        power -= ramiPower;
    }

    generateLighting2d(verts, inds, params, midP, dst, power, level + 1, midP_ind);
}

static const vec4 k_transparent = { 1, 1, 1, 0 };
static const vec4 k_opaque = { 1, 1, 1, 1 };

static void lightning_linesToTriangles_rec(
    Lightning_linesToTriangles_outParams& out,
    std::vector<u32>& toReorient,
    const Lightning_linesToTriangles_inParams& in,
    const std::vector<Childs>& m,
    u32 parent, float parentR, u32 node, u32 parentI0, u32 parentI1, u32 parentI2)
{
    auto to3d = [h = in.params.height](vec2 p)
    {
        return glm::vec3(p.x, h, -p.y);
    };
    auto doQuad = [&out](u32 i0, u32 i1, u32 i2, u32 i3)
    {
        out.inds.push_back(i0);
        out.inds.push_back(i1);
        out.inds.push_back(i2);
        out.inds.push_back(i0);
        out.inds.push_back(i2);
        out.inds.push_back(i3);
    };

    const size_t numChildren = m[node].size();
    if (numChildren == 0) {
        const u32 i0 = out.verts.size();
        out.verts.push_back({ to3d(in.verts[node].pos), k_transparent });
        out.inds.push_back(i0);
        out.inds.push_back(parentI1);
        out.inds.push_back(parentI0);
        out.inds.push_back(i0);
        out.inds.push_back(parentI2);
        out.inds.push_back(parentI1);
    }
    else {
        const float& power = in.verts[node].power;
        const float r = in.params.maxThickness * sqrt(power);
        const vec2& a = in.verts[parent].pos;
        const vec2& b = in.verts[node].pos;
        vec2 ab = b - a;
        const float abLen = length(ab);
        ab /= abLen;
        if (numChildren == 1) {
            const u32 child = m[node][0];
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

            if (!isnan(in.params.cameraPos.x)) {
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
            assert(numChildren == 2);

            u32 childL = m[node][0];
            u32 childR = m[node][1];

            const vec2 abRight = { -ab.y, b.x };
            if (dot(abRight, normalize(in.verts[childL].pos - b)) > dot(abRight, normalize(in.verts[childR].pos - b)))
                std::swap(childL, childR);
            const vec2& cl = in.verts[childL].pos;
            const vec2& cr = in.verts[childR].pos;
            const float rl = in.params.maxThickness * sqrt(in.verts[childL].power);
            const float rr = in.params.maxThickness * sqrt(in.verts[childR].power);

            vec2 bcl = cl - b;
            const float bclLen = length(bcl);
            bcl /= bclLen;
            vec2 bcr = cr - b;
            const float bcrLen = length(bcr);
            bcr /= bcrLen;

            const vec2 bclRight = { -bcl.y, bcl.x };
            const vec2 bcrRight = { -bcr.y, bcr.x };

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
                    d = glm::clamp(d, 0.f, 1.2f * abLen);
                    leftP = parentL + d * ab;
                }
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
                    d = glm::clamp(d, 0.f, 1.2f * abLen);
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
            }

            u32 i0 = out.verts.size();
            out.verts.push_back({ to3d(leftP), k_transparent });
            out.verts.push_back({ to3d(b), k_opaque });
            out.verts.push_back({ to3d(rightP), k_transparent });
            out.verts.push_back({ to3d(midP), k_transparent });

            if (!isnan(in.params.cameraPos.x)) {
                toReorient.push_back(i0 + 1);
                toReorient.push_back(i0);
                toReorient.push_back(i0 + 1);
                toReorient.push_back(i0 + 2);
                toReorient.push_back(i0 + 1);
                toReorient.push_back(i0 + 3);
            }

            doQuad(i0, i0 + 1, parentI1, parentI0);
            doQuad(i0 + 1, i0 + 2, parentI2, parentI1);

            lightning_linesToTriangles_rec(out, toReorient, in, m, node, rl, childL, i0, i0 + 1, i0 + 3);
            lightning_linesToTriangles_rec(out, toReorient, in, m, node, rr, childR, i0 + 3, i0 + 1, i0 + 2);
        }
    }
}

void lightning_linesToTriangles(
    Lightning_linesToTriangles_outParams& out,
    const Lightning_linesToTriangles_inParams& in)
{
    auto to3d = [h = in.params.height](vec2 p) {
        return vec3(p.x, h, -p.y);
    };
    std::vector<Childs> m(in.verts.size());
    for (size_t i = 0; i < in.inds.size(); i += 2) {
        m[in.inds[i]].push_back(in.inds[i + 1]);
    }

    std::vector<bool> isRoot(in.verts.size(), true);
    for (u32 ori = 0; ori < isRoot.size(); ori++) {
        for (u32 i = 0; i < 2; i++) {
            u32 dst = m[ori][i];
            if (dst == -1)
                break;
            isRoot[dst] = false;
        }
    }

    for(u32 root = 0; root < isRoot.size(); root++)
    {
        if (!isRoot[root])
            continue;
        assert(m[root].size() == 1);

        const u32 rootChild = m[root][0];
        const vec2& a = in.verts[root].pos;
        const vec2& b = in.verts[rootChild].pos;
        const vec2 ab = normalize(b - a);
        const vec2 right = normalize(vec2(-ab.y, ab.x));
        const float& power = in.verts[root].power;
        const float r = in.params.maxThickness * sqrt(power);
        const vec2 al = a - right * r;
        const vec2 ar = a + right * r;
        const u32 i0 = out.verts.size();
        out.verts.push_back({ to3d(al), k_transparent });
        out.verts.push_back({ to3d(a), k_opaque });
        out.verts.push_back({ to3d(ar), k_transparent });
        std::vector<u32> toReorient;
        if (!isnan(in.params.cameraPos.x)) {
            toReorient.reserve(in.inds.size() * 2);
            toReorient.push_back(i0 + 1);
            toReorient.push_back(i0 + 0);
            toReorient.push_back(i0 + 1);
            toReorient.push_back(i0 + 2);
        }
        lightning_linesToTriangles_rec(out, toReorient, in, m, root, r,
            rootChild, i0, i0 + 1, i0 + 2);

        for (size_t i = 0; i < toReorient.size(); i += 2) {
            const vec3& midP = out.verts[toReorient[i]].pos;
            vec3& p = out.verts[toReorient[i + 1]].pos;

            const vec3 toCam = normalize(in.params.cameraPos - midP);
            vec3 toP = p - midP;
            const float r = length(toP);
            toP = toP - toCam * dot(toCam, toP);
            toP = r * normalize(toP);
            p = midP + toP;
        }
    }
}

void generateLightning_triangles(
    std::vector<LightningVert3d>& verts, std::vector<u32>& inds,
    const LightningParams& params,
    const glm::vec2& ori, const glm::vec2& dst
)
{
    std::vector<LightningVert2d> lines_verts;
    std::vector<u32> lines_inds;
    generateLighting2d(lines_verts, lines_inds, params, ori, dst);

    Lightning_linesToTriangles_outParams outParams{ verts, inds };
    const Lightning_linesToTriangles_inParams inParams{
        lines_verts,
        lines_inds,
        params
    };

    lightning_linesToTriangles(outParams, inParams);
}

void generateLightning_lines(
    std::vector<LightningVert3d>& verts, std::vector<u32>& inds,
    const LightningParams& params,
    const glm::vec2& ori, const glm::vec2& dst
)
{
    std::vector<LightningVert2d> verts2d;
    generateLighting2d(verts2d, inds, params, ori, dst);

    lightning_linesTo3d(verts, verts2d);
}

void lightning_linesTo3d(
    std::vector<LightningVert3d>& out,
    std::span<const LightningVert2d> in)
{
    out.reserve(out.size() + in.size());
    for (auto& v : in) {
        out.push_back({ to3d(v.pos), {1, 1, 0, 1} });
    }
}