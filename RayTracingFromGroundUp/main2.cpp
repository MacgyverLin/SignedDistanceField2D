#include "bitmap.h"
#include "vec2.h"
#include "vec3.h"
#include "vec4.h"
#include "mat2.h"

float sdCircle(vec2 p, float r)
{
	return length(p) - r;
}

float sdRoundedBox(vec2 p, vec2 b, const vec4& radius)
{
	vec4 r = radius;

	r.x = (p.x > 0.0) ? r.x : r.w;
	r.y = (p.x > 0.0) ? r.y : r.z;

	r.x = (p.y > 0.0) ? r.x : r.y;
	vec2 q = abs(p) - b + vec2(r.x);
	return min(max(q.x, q.y), 0.0f) + length(max(q, 0.0f)) - r.x;
}

float sdBox(vec2 p, vec2 b)
{
	vec2 d = abs(p) - b;
	
	return length(max(d, 0.0f)) + min(max(d.x, d.y), 0.0f);
}

float sdOrientedBox(vec2 p, vec2 a, vec2 b, float th)
{
	float l = length(b - a);
	vec2 d = (b - a) / l;
	vec2 q = (p - (a + b) * 0.5);
	q = mat2(d.x, -d.y, d.y, d.x) * q;
	q = abs(q) - vec2(l, th) * 0.5;
	return length(max(q, 0.0)) + min(max(q.x, q.y), 0.0f);
}

float sdSegment(vec2 p, vec2 a, vec2 b)
{
	vec2 pa = p - a;
	vec2 ba = b - a;
	float h = clamp(dot(pa, ba) / dot(ba, ba), 0.0f, 1.0f);
	return length(pa - ba * h);
}

float sdRhombus(vec2 p, vec2 b)
{
	vec2 q = abs(p);
	float h = clamp((-2.0f * ndot(q, b) + ndot(b, b)) / dot(b, b), -1.0f, 1.0f);
	float d = length(q - 0.5 * b * vec2(1.0 - h, 1.0 + h));
	return d * sign(q.x * b.y + q.y * b.x - b.x * b.y);
}

float sdEquilateralTriangle(vec2 p)
{
	const float k = sqrt(3.0f);
	p.x = abs(p.x) - 1.0f;
	p.y = p.y + 1.0f / k;
	if (p.x + k * p.y > 0.0f) 
		p = vec2(p.x - k * p.y, -k * p.x - p.y) / 2.0f;
	p.x -= clamp(p.x, -2.0f, 0.0f);
	return -length(p) * sign(p.y);
}

float sdTriangleIsosceles(vec2 p, vec2 q)
{
	p.x = abs(p.x);
	vec2 a = p - q * clamp(dot(p, q) / dot(q, q), 0.0f, 1.0f);
	vec2 b = p - q * vec2(clamp(p.x / q.x, 0.0f, 1.0f), 1.0f);
	float s = -sign(q.y);
	vec2 d = min(vec2(dot(a, a), s * (p.x * q.y - p.y * q.x)), vec2(dot(b, b), s * (p.y - q.y)));
	return -sqrt(d.x) * sign(d.y);
}

float sdUnevenCapsule(vec2 p, float r1, float r2, float h)
{
	p.x = abs(p.x);
	float b = (r1 - r2) / h;
	float a = sqrt(1.0f - b * b);
	float k = dot(p, vec2(-b, a));
	if (k < 0.0) return length(p) - r1;
	if (k > a * h) return length(p - vec2(0.0f, h)) - r2;
	return dot(p, vec2(a, b)) - r1;
}

float sdTriangle(vec2 p, vec2 p0, vec2 p1, vec2 p2)
{
	vec2 e0 = p1 - p0, e1 = p2 - p1, e2 = p0 - p2;
	vec2 v0 = p - p0, v1 = p - p1, v2 = p - p2;
	vec2 pq0 = v0 - e0 * clamp(dot(v0, e0) / dot(e0, e0), 0.0f, 1.0f);
	vec2 pq1 = v1 - e1 * clamp(dot(v1, e1) / dot(e1, e1), 0.0f, 1.0f);
	vec2 pq2 = v2 - e2 * clamp(dot(v2, e2) / dot(e2, e2), 0.0f, 1.0f);
	float s = sign(e0.x * e2.y - e0.y * e2.x);
	vec2 d = min(min(vec2(dot(pq0, pq0), s * (v0.x * e0.y - v0.y * e0.x)),
		vec2(dot(pq1, pq1), s * (v1.x * e1.y - v1.y * e1.x))),
		vec2(dot(pq2, pq2), s * (v2.x * e2.y - v2.y * e2.x)));
	return -sqrt(d.x) * sign(d.y);
}

float sdPentagon(vec2 p, float r)
{
    const vec3 k = vec3(0.809016994, 0.587785252, 0.726542528);
    p.x = abs(p.x);
    p -= 2.0 * min(dot(vec2(-k.x, k.y), p), 0.0f) * vec2(-k.x, k.y);
    p -= 2.0 * min(dot(vec2(k.x, k.y), p), 0.0f) * vec2(k.x, k.y);
    p -= vec2(clamp(p.x, -r * k.z, r * k.z), r);
    return length(p) * sign(p.y);
}

float sdHexagon(vec2 p, float r)
{
    const vec3 k = vec3(-0.866025404, 0.5, 0.577350269);
    p = abs(p);
    p -= 2.0 * min(dot(vec2(k.x, k.y), p), 0.0f) * vec2(k.x, k.y);
    p -= vec2(clamp(p.x, -k.z * r, k.z * r), r);
    return length(p) * sign(p.y);
}

float sdOctogon(vec2 p, float r)
{
    const vec3 k = vec3(-0.9238795325, 0.3826834323, 0.4142135623);
    p = abs(p);
    p -= 2.0 * min(dot(vec2(k.x, k.y), p), 0.0f) * vec2(k.x, k.y);
    p -= 2.0 * min(dot(vec2(-k.x, k.y), p), 0.0f) * vec2(-k.x, k.y);
    p -= vec2(clamp(p.x, -k.z * r, k.z * r), r);
    return length(p) * sign(p.y);
}

float sdHexagram(vec2 p, float r)
{
    const vec4 k = vec4(-0.5, 0.8660254038, 0.5773502692, 1.7320508076);
    p = abs(p);
    p -= 2.0 * min(dot(vec2(k.x, k.y), p), 0.0f) * vec2(k.x, k.y);
    p -= 2.0 * min(dot(vec2(k.y, k.x), p), 0.0f) * vec2(k.y, k.x);
    p -= vec2(clamp(p.x, r * k.z, r * k.w), r);
    return length(p) * sign(p.y);
}

float sdStar5(vec2 p, float r, float rf)
{
    const vec2 k1 = vec2(0.809016994375, -0.587785252292);
    const vec2 k2 = vec2(-k1.x, k1.y);
    p.x = abs(p.x);
    p -= 2.0 * max(dot(k1, p), 0.0f) * k1;
    p -= 2.0 * max(dot(k2, p), 0.0f) * k2;
    p.x = abs(p.x);
    p.y -= r;
    vec2 ba = rf * vec2(-k1.y, k1.x) - vec2(0, 1);
    float h = clamp(dot(p, ba) / dot(ba, ba), 0.0f, r);
    return length(p - ba * h) * sign(p.y * ba.x - p.x * ba.y);
}

float sdStar(vec2 p, float r, int n, float m)
{
    // these 4 lines can be precomputed for a given shape
    float an = 3.141593 / float(n);
    float en = 3.141593 / m;
    vec2  acs = vec2(cos(an), sin(an));
    vec2  ecs = vec2(cos(en), sin(en)); // ecs=vec2(0,1) and simplify, for regular polygon,

    // reduce to first sector
    float bn = fmod(atan2(p.x, p.y), 2.0f * an) - an;
    p = length(p) * vec2(cos(bn), abs(sin(bn)));

    // line sdf
    p -= r * acs;
    p += ecs * clamp(-dot(p, ecs), 0.0f, r * acs.y / ecs.y);
    return length(p) * sign(p.x);
}

float dot2(vec2 v) { return dot(v, v); }

float sdTrapezoid(vec2 p, float r1, float r2, float he)
{
    vec2 k1 = vec2(r2, he);
    vec2 k2 = vec2(r2 - r1, 2.0 * he);
    p.x = abs(p.x);
    vec2 ca = vec2(p.x - min(p.x, (p.y < 0.0) ? r1 : r2), abs(p.y) - he);
    vec2 cb = p - k1 + k2 * clamp(dot(k1 - p, k2) / dot2(k2), 0.0f, 1.0f);
    float s = (cb.x < 0.0 && ca.y < 0.0) ? -1.0 : 1.0;
    return s * sqrt(min(dot2(ca), dot2(cb)));
}

float sdPie(vec2 p, vec2 c, float r)
{
    p.x = abs(p.x);
    float l = length(p) - r;
    float m = length(p - c * clamp(dot(p, c), 0.0f, r)); // c = sin/cos of the aperture
    return max(l, m * sign(c.y * p.x - c.x * p.y));
}

float sdArc(vec2 p, vec2 sca, vec2 scb, float ra, float rb)
{
    p = mat2(sca.x, sca.y, -sca.y, sca.x) * p;
    p.x = abs(p.x);
    float k = (scb.y * p.x > scb.x * p.y) ? dot(vec2(p.x, p.y), scb) : length(vec2(p.x, p.y));
    return sqrt(dot(p, p) + ra * ra - 2.0 * ra * k) - rb;
}

float sdHorseshoe(vec2 p, vec2 c, float r, vec2 w)
{
    p.x = abs(p.x);
    float l = length(p);
    p = mat2(-c.x, c.y,
        c.y, c.x) * p;
    p = vec2((p.y > 0.0) ? p.x : l * sign(-c.x),
        (p.x > 0.0) ? p.y : l);
    p = vec2(p.x, abs(p.y - r)) - w;
    return length(max(p, 0.0)) + min(0.0f, max(p.x, p.y));
}

float sdVesica(vec2 p, float r, float d)
{
    p = abs(p);
    float b = sqrt(r * r - d * d);
    return ((p.y - b) * d > p.x * b) ? length(p - vec2(0.0, b)) : length(p - vec2(-d, 0.0)) - r;
}

float sdEgg(vec2 p, float ra, float rb)
{
    const float k = sqrt(3.0);
    p.x = abs(p.x);
    float r = ra - rb;
    return ((p.y < 0.0) ? length(vec2(p.x, p.y)) - r :
        (k * (p.x + r) < p.y) ? length(vec2(p.x, p.y - k * r)) :
        length(vec2(p.x + r, p.y)) - 2.0 * r) - rb;
}

float sdCross(vec2 p, vec2 b, float r)
{
    p = abs(p); p = (p.y > p.x) ? vec2(p.y, p.x) : vec2(p.x, p.y);
    vec2  q = p - b;
    float k = max(q.y, q.x);
    vec2  w = (k > 0.0) ? q : vec2(b.y - p.x, -k);
    return sign(k) * length(max(w, 0.0)) + r;
}

float sdRoundedX(vec2 p, float w, float r)
{
    p = abs(p);
    return length(p - min(p.x + p.y, w) * 0.5) - r;
}

/*
float sdPolygon(vec2* v, int N, vec2 p)
{
    float d = dot(p - v[0], p - v[0]);
    float s = 1.0;
    for (int i = 0, j = N - 1; i < N; j = i, i++)
    {
        vec2 e = v[j] - v[i];
        vec2 w = p - v[i];
        vec2 b = w - e * clamp(dot(w, e) / dot(e, e), 0.0f, 1.0f);
        d = min(d, dot(b, b));
        bvec3 c = bvec3(p.y >= v[i].y, p.y<v[j].y, e.x* w.y>e.y * w.x);
        if (all(c) || all(not(c)))
            s *= -1.0;
    }
    return s * sqrt(d);
}
*/

float sdEllipse(vec2 p, vec2 ab)
{
    p = abs(p); 
    if (p.x > p.y) 
    { 
        p = vec2(p.y, p.x); 
        ab = vec2(ab.y, ab.x); 
    }
    float l = ab.y * ab.y - ab.x * ab.x;
    float m = ab.x * p.x / l;      float m2 = m * m;
    float n = ab.y * p.y / l;      float n2 = n * n;
    float c = (m2 + n2 - 1.0) / 3.0; float c3 = c * c * c;
    float q = c3 + m2 * n2 * 2.0;
    float d = c3 + m2 * n2;
    float g = m + m * n2;
    float co;
    if (d < 0.0)
    {
        float h = acos(q / c3) / 3.0;
        float s = cos(h);
        float t = sin(h) * sqrt(3.0);
        float rx = sqrt(-c * (s + t + 2.0) + m2);
        float ry = sqrt(-c * (s - t + 2.0) + m2);
        co = (ry + sign(l) * rx + abs(g) / (rx * ry) - m) / 2.0;
    }
    else
    {
        float h = 2.0 * m * n * sqrt(d);
        float s = sign(q + h) * pow(abs(q + h), 1.0 / 3.0);
        float u = sign(q - h) * pow(abs(q - h), 1.0 / 3.0);
        float rx = -s - u - c * 4.0 + 2.0 * m2;
        float ry = (s - u) * sqrt(3.0);
        float rm = sqrt(rx * rx + ry * ry);
        co = (ry / sqrt(rm - rx) + 2.0 * g / rm - m) / 2.0;
    }
    vec2 r = ab * vec2(co, sqrt(1.0 - co * co));
    return length(r - p) * sign(p.y - r.y);
}

float sdParabola(vec2 pos, float k)
{
    pos.x = abs(pos.x);
    float ik = 1.0 / k;
    float p = ik * (pos.y - 0.5 * ik) / 3.0;
    float q = 0.25 * ik * ik * pos.x;
    float h = q * q - p * p * p;
    float r = sqrt(abs(h));
    float x = (h > 0.0) ?
        pow(q + r, 1.0 / 3.0) - pow(abs(q - r), 1.0 / 3.0) * sign(r - q) :
        2.0 * cos(atan2(r, q) / 3.0) * sqrt(p);
    return length(pos - vec2(x, k * x * x)) * sign(pos.x - x);
}

float sdParabola(vec2 pos, float wi, float he)
{
    pos.x = abs(pos.x);
    float ik = wi * wi / he;
    float p = ik * (he - pos.y - 0.5 * ik) / 3.0;
    float q = pos.x * ik * ik * 0.25;
    float h = q * q - p * p * p;
    float r = sqrt(abs(h));
    float x = (h > 0.0) ?
        pow(q + r, 1.0 / 3.0) - pow(abs(q - r), 1.0 / 3.0) * sign(r - q) :
        2.0 * cos(atan(r / q) / 3.0) * sqrt(p);
    x = min(x, wi);
    return length(pos - vec2(x, he - x * x / ik)) *
        sign(ik * (pos.y - he) + pos.x * pos.x);
}

float sdBezier(vec2 pos, vec2 A, vec2 B, vec2 C)
{
    vec2 a = B - A;
    vec2 b = A - 2.0 * B + C;
    vec2 c = a * 2.0;
    vec2 d = A - pos;
    float kk = 1.0 / dot(b, b);
    float kx = kk * dot(a, b);
    float ky = kk * (2.0 * dot(a, a) + dot(d, b)) / 3.0;
    float kz = kk * dot(d, a);
    float res = 0.0;
    float p = ky - kx * kx;
    float p3 = p * p * p;
    float q = kx * (2.0 * kx * kx - 3.0 * ky) + kz;
    float h = q * q + 4.0 * p3;
    if (h >= 0.0)
    {
        h = sqrt(h);
        vec2 x = (vec2(h, -h) - q) / 2.0;
        vec2 uv = sign(x) * pow(abs(x), vec2(1.0f / 3.0f));
        float t = clamp(uv.x + uv.y - kx, 0.0f, 1.0f);
        res = dot2(d + (c + b * t) * t);
    }
    else
    {
        float z = sqrt(-p);
        float v = acos(q / (p * z * 2.0)) / 3.0;
        float m = cos(v);
        float n = sin(v) * 1.732050808;
        vec3  t = clamp(vec3(m + m, -n - m, n - m) * z - kx, 0.0, 1.0);
        res = min(dot2(d + (c + b * t.x) * t.x),
            dot2(d + (c + b * t.y) * t.y));
        // the third root cannot be the closest
        // res = min(res,dot2(d+(c+b*t.z)*t.z));
    }
    return sqrt(res);
}

float opRound(vec2 p, std::function<float(vec2)> sdFunc, float r)
{
    return sdFunc(p) - r;
}

float opOnion(vec2 p, std::function<float(vec2)> sdFunc, float r)
{
    return abs(sdFunc(p)) - r;
}

float opUnion(vec2 p, std::function<float(vec2)> d1, std::function<float(vec2)> d2)
{ 
    return min(d1(p), d2(p));
}

float opSubtraction(vec2 p, std::function<float(vec2)> d1, std::function<float(vec2)> d2)
{ 
    return max(-d1(p), d2(p));
}

float opIntersection(vec2 p, std::function<float(vec2)> d1, std::function<float(vec2)> d2)
{ 
    return max(d1(p), d2(p)); 
}

float opSmoothUnion(vec2 p, std::function<float(vec2)> d1, std::function<float(vec2)> d2, float k) 
{
    float h = clamp(0.5 + 0.5 * (d2(p) - d1(p)) / k, 0.0, 1.0);
    return mix(d2(p), d1(p), h) - k * h * (1.0 - h);
}

float opSmoothSubtraction(vec2 p, std::function<float(vec2)> d1, std::function<float(vec2)> d2, float k) 
{
    float h = clamp(0.5 - 0.5 * (d2(p) + d1(p)) / k, 0.0, 1.0);
    return mix(d2(p), -d1(p), h) + k * h * (1.0 - h);
}

float opSmoothIntersection(vec2 p, std::function<float(vec2)> d1, std::function<float(vec2)> d2, float k) 
{
    float h = clamp(0.5 - 0.5 * (d2(p) - d1(p)) / k, 0.0, 1.0);
    return mix(d2(p), d1(p), h) + k * h * (1.0 - h);
}

vec3 color(float d)
{
	vec3 col = vec3(1.0f) - sign(d) * vec3(0.1, 0.4, 0.7);
	col *= (1.0 - exp(-2.0 * abs(d)));
	col *= (0.8 + 0.2 * cos(140.0 * d));
	col = mix(col, vec3(1.0), 1.0 - smoothstep(0.0f, 0.02f, abs(d)));

	return col;
}

void draw(int w, int h, int N, const char* filename, std::function<float(vec2 p)> sdfFunc)
{
	bitmap bmp(w, h);

	int sqrtN = sqrt(N);
	float ddu = 1.0f / sqrtN;
	float ddv = 1.0f / sqrtN;

	for (int y = 0; y < h; y++)
	{
		for (int x = 0; x < w; x++)
		{
			vec3 col(0, 0, 0);
			vec2 p = (vec2(x, y) * 2.0 - vec2(w, h)) / h;
			// printf("%f, %f\n", p.x, p.y);
			float distance = sdfFunc(vec2(p.x, p.y));
			col = color(distance);

			bmp.setPixel(x, y, col);
		}
	}

	bmp.tonemap();
	bmp.saveBMP(filename);
}

float opTx(vec2 p, mat2 t, std::function<float(vec2)> sdfFunc)
{
    return sdfFunc(inverse(t) * p);
}

float opScale(vec2 p, float s, std::function<float(vec2)> sdfFunc)
{
    return sdfFunc(p / s) * s;
}

float opSymX(vec2 p, std::function<float(vec2)> sdfFunc)
{
    p.x = abs(p.x);
    return sdfFunc(p);
}

float opSymY(vec2 p, std::function<float(vec2)> sdfFunc)
{
    p.y = abs(p.y);
    return sdfFunc(p);
}

float opRep(vec2 p, vec2 c, std::function<float(vec2)> sdfFunc)
{
    vec2 q = mod(p + 0.5 * c, c) - 0.5 * c;
    return sdfFunc(q);
}

int main()
{
	int image_width = 256;
	int image_height = 256;
    
    float radius = 0.0f;

    std::function<float(vec2 p)> sdf = nullptr;
    std::function<float(vec2 p)> sdf1 = nullptr;
    std::function<float(vec2 p)> sdf2 = nullptr;
    std::function<float(vec2 p)> modifier = nullptr;

    sdf1 = bind(sdCircle, std::placeholders::_1, 0.6f);
    sdf2 = bind(sdTriangleIsosceles, std::placeholders::_1, vec2(0.4, 0.7));
    sdf = bind(opUnion, std::placeholders::_1, sdf1, sdf2);
    modifier = bind(opOnion, std::placeholders::_1, sdf, radius);
    draw(image_width, image_height, 1, "opUnion.bmp", sdf);

    sdf1 = bind(sdCircle, std::placeholders::_1, 0.6f);
    sdf2 = bind(sdTriangleIsosceles, std::placeholders::_1, vec2(0.4, 0.7));
    sdf = bind
    (
        opSubtraction, 
        std::placeholders::_1, 
        sdf1,
        sdf2
    );
    modifier = bind(opOnion, std::placeholders::_1, sdf, radius);
    draw(image_width, image_height, 1, "opSubtractions.bmp", sdf);

    sdf1 = bind(sdCircle, std::placeholders::_1, 0.6f);
    sdf2 = bind(sdTriangleIsosceles, std::placeholders::_1, vec2(0.4, 0.7));
    sdf = bind(opIntersection, std::placeholders::_1, sdf1, sdf2);
    modifier = bind(opOnion, std::placeholders::_1, sdf, radius);
    draw(image_width, image_height, 1, "opIntersection.bmp", sdf);

    mat2 rot(0.0f, 0.0f, 0.0f, 0.0f);
    rot.make_rotate(30);
    sdf1 = bind(sdUnevenCapsule, std::placeholders::_1, 0.2, 0.3, 0.5);
    sdf = bind(opTx, std::placeholders::_1, rot, sdf1);
    draw(image_width, image_height, 1, "opTx.bmp", sdf);

    sdf1 = bind(sdUnevenCapsule, std::placeholders::_1, 0.2, 0.3, 0.5);
    sdf = bind(opScale, std::placeholders::_1, 2.0f, sdf1);
    draw(image_width, image_height, 1, "opScale.bmp", sdf);

    sdf1 = bind(sdUnevenCapsule, std::placeholders::_1, 0.2, 0.3, 0.5);
    sdf = bind(opSymX, std::placeholders::_1, sdf1);
    draw(image_width, image_height, 1, "opSymX.bmp", sdf);

    sdf1 = bind(sdUnevenCapsule, std::placeholders::_1, 0.2, 0.3, 0.5);
    sdf = bind(opSymY, std::placeholders::_1, sdf1);
    draw(image_width, image_height, 1, "opSymY.bmp", sdf);

    sdf1 = bind(sdBox, std::placeholders::_1, vec2(0.02, 0.02));
    sdf = bind(opRep, std::placeholders::_1, vec2(0.2), sdf1);
    draw(image_width, image_height, 1, "opRep.bmp", sdf);    


    sdf1 = bind(sdCircle, std::placeholders::_1, 0.6f);
    sdf2 = bind(sdTriangleIsosceles, std::placeholders::_1, vec2(0.4, 0.7));
    sdf = bind(opIntersection, std::placeholders::_1, sdf1, sdf2);
    modifier = bind(opOnion, std::placeholders::_1, sdf, radius);
    draw(image_width, image_height, 1, "opIntersection.bmp", sdf);

    sdf1 = bind(sdCircle, std::placeholders::_1, 0.6f);
    sdf2 = bind(sdTriangleIsosceles, std::placeholders::_1, vec2(0.4, 0.7));
    sdf = bind(opIntersection, std::placeholders::_1, sdf1, sdf2);
    modifier = bind(opOnion, std::placeholders::_1, sdf, radius);
    draw(image_width, image_height, 1, "opIntersection.bmp", sdf);

    sdf = bind(sdCircle, std::placeholders::_1, 0.3f);
    modifier = bind(opOnion, std::placeholders::_1, sdf, radius);
    draw(image_width, image_height, 1, "sdCircle.bmp", modifier);

    sdf = bind(sdRoundedBox, std::placeholders::_1, vec2(0.5, 0.5), vec4(0.03, 0.06, 0.09, 0.12));
    modifier = bind(opOnion, std::placeholders::_1, sdf, radius);
	draw(image_width, image_height, 1, "sdRoundedBox.bmp", bind(opOnion, std::placeholders::_1, modifier, radius));
	
    sdf = bind(sdBox, std::placeholders::_1, vec2(0.5, 0.5));
    modifier = bind(opOnion, std::placeholders::_1, sdf, radius);
    draw(image_width, image_height, 1, "sdBox.bmp", bind(opOnion, std::placeholders::_1, modifier, radius) );
	    
    sdf = bind(sdTriangleIsosceles, std::placeholders::_1, vec2(0.4, 0.7));
    modifier = bind(opOnion, std::placeholders::_1, sdf, radius);
    draw(image_width, image_height, 1, "sdTriangleIsosceles.bmp", bind(opOnion, std::placeholders::_1, modifier, radius));

    sdf = bind(sdOrientedBox, std::placeholders::_1, vec2(0.03, 0.07), vec2(0.09, 0.05), 0);
    modifier = bind(opOnion, std::placeholders::_1, sdf, radius);
    draw(image_width, image_height, 1, "sdOrientedBox.bmp", bind(opOnion, std::placeholders::_1, modifier, radius));
	
    sdf = bind(sdSegment, std::placeholders::_1, vec2(-0.3, -0.3), vec2(0.3, 0.3));
    modifier = bind(opOnion, std::placeholders::_1, sdf, radius);
    draw(image_width, image_height, 1, "sdSegment.bmp", bind(opOnion, std::placeholders::_1, modifier, radius));
	
    sdf = bind(sdRhombus, std::placeholders::_1, vec2(0.4, 0.7));
    modifier = bind(opOnion, std::placeholders::_1, sdf, radius);
    draw(image_width, image_height, 1, "sdRhombus.bmp", bind(opOnion, std::placeholders::_1, modifier, radius));
	
    sdf = bind(sdEquilateralTriangle, std::placeholders::_1);
    modifier = bind(opOnion, std::placeholders::_1, sdf, radius);
    draw(image_width, image_height, 1, "sdEquilateralTriangle.bmp", bind(opOnion, std::placeholders::_1, modifier, radius));
	
    sdf = bind(sdTriangleIsosceles, std::placeholders::_1, vec2(0.4, 0.7));
    modifier = bind(opOnion, std::placeholders::_1, sdf, radius);
    draw(image_width, image_height, 1, "sdTriangleIsosceles.bmp", bind(opOnion, std::placeholders::_1, modifier, radius));
	
    sdf = bind(sdUnevenCapsule, std::placeholders::_1, 0.2, 0.3, 0.5);
    modifier = bind(opOnion, std::placeholders::_1, sdf, radius);
    draw(image_width, image_height, 1, "sdUnevenCapsule.bmp", bind(opOnion, std::placeholders::_1, modifier, radius));
    
	vec2 v1 = cos(vec2(0.0, 2.00) + 0.0);
	vec2 v2 = cos(vec2(0.0, 1.50) + 1.5);
	vec2 v3 = cos(vec2(0.0, 3.00) + 4.0);
    sdf = bind(sdTriangle, std::placeholders::_1, v1, v2, v3);
    modifier = bind(opOnion, std::placeholders::_1, sdf, radius);
	draw(image_width, image_height, 1, "sdTriangle.bmp", bind(opOnion, std::placeholders::_1, modifier, radius));
    
    sdf = bind(sdPentagon, std::placeholders::_1, 0.5);
    modifier = bind(opOnion, std::placeholders::_1, sdf, radius);
    draw(image_width, image_height, 1, "sdPentagon.bmp", bind(opOnion, std::placeholders::_1, modifier, radius));
    
    sdf = bind(sdHexagon, std::placeholders::_1, 0.5);
    modifier = bind(opOnion, std::placeholders::_1, sdf, radius);
    draw(image_width, image_height, 1, "sdHexagon.bmp", bind(opOnion, std::placeholders::_1, modifier, radius));
    
    sdf = bind(sdOctogon, std::placeholders::_1, 0.5);
    modifier = bind(opOnion, std::placeholders::_1, sdf, radius);
    draw(image_width, image_height, 1, "sdOctogon.bmp", bind(opOnion, std::placeholders::_1, modifier, radius));
    
    sdf = bind(sdHexagram, std::placeholders::_1, 0.5);
    modifier = bind(opOnion, std::placeholders::_1, sdf, radius);
    draw(image_width, image_height, 1, "sdHexagram.bmp", bind(opOnion, std::placeholders::_1, modifier, radius));
    
    sdf = bind(sdStar5, std::placeholders::_1, 0.5, 0.8);
    modifier = bind(opOnion, std::placeholders::_1, sdf, radius);
    draw(image_width, image_height, 1, "sdStar5.bmp", bind(opOnion, std::placeholders::_1, modifier, radius));

    //float t = 10.0f / 3.0;
    //float n = 3.0 + fmod(floor(t), 9.0);      // n, number of sides
    //float a = t - floor(t);                   // angle factor
    //float m = 2.0f + a * a * (n - 2.0f);       // angle divisor, between 2 and n    
    sdf = bind(sdStar5, std::placeholders::_1, 0.5, 0.8);
    modifier = bind(opOnion, std::placeholders::_1, sdf, radius);
    draw(image_width, image_height, 1, "sdStar.bmp", bind(opOnion, std::placeholders::_1, modifier, radius));
    
    sdf = bind(sdTrapezoid, std::placeholders::_1, 0.5, 0.8, 1);
    modifier = bind(opOnion, std::placeholders::_1, sdf, radius);
    draw(image_width, image_height, 1, "sdTrapezoid.bmp", bind(opOnion, std::placeholders::_1, modifier, radius));

    sdf = bind(sdPie, std::placeholders::_1, vec2(0, 0), 0.5);
    modifier = bind(opOnion, std::placeholders::_1, sdf, radius);
    draw(image_width, image_height, 1, "sdPie.bmp", bind(opOnion, std::placeholders::_1, modifier, radius));

    sdf = bind(sdArc, std::placeholders::_1, vec2(-1, 0), vec2(1, 0), 0.5, 0.7);
    modifier = bind(opOnion, std::placeholders::_1, sdf, radius);
    draw(image_width, image_height, 1, "sdArc.bmp", bind(opOnion, std::placeholders::_1, modifier, radius));

    sdf = bind(sdHorseshoe, std::placeholders::_1, vec2(0, 0), 0.5, vec2(0.5, 0.3));
    modifier = bind(opOnion, std::placeholders::_1, sdf, radius);
    draw(image_width, image_height, 1, "sdHorseshoe.bmp", bind(opOnion, std::placeholders::_1, modifier, radius));

    sdf = bind(sdVesica, std::placeholders::_1, 0.2, 0.7);
    modifier = bind(opOnion, std::placeholders::_1, sdf, radius);
    draw(image_width, image_height, 1, "sdVesica.bmp", bind(opOnion, std::placeholders::_1, modifier, radius));

    sdf = bind(sdEgg, std::placeholders::_1, 0.2, 0.7);
    modifier = bind(opOnion, std::placeholders::_1, sdf, radius);
    draw(image_width, image_height, 1, "sdEgg.bmp", bind(opOnion, std::placeholders::_1, modifier, radius));

    sdf = bind(sdCross, std::placeholders::_1, vec2(0.5, 0.5), 0.7);
    modifier = bind(opOnion, std::placeholders::_1, sdf, radius);
    draw(image_width, image_height, 1, "sdCross.bmp", bind(opOnion, std::placeholders::_1, modifier, radius));

    sdf = bind(sdRoundedX, std::placeholders::_1, 0.6, 0.6);
    modifier = bind(opOnion, std::placeholders::_1, sdf, radius);
    draw(image_width, image_height, 1, "sdRoundedX.bmp", bind(opOnion, std::placeholders::_1, modifier, radius));

	return -1;
}

/*
bind(sdPolygon(vec2* v, int N, std::placeholders::_1)
{
    float d = dot(p - v[0], p - v[0]);
    float s = 1.0;
    for (int i = 0, j = N - 1; i < N; j = i, i++)
    {
        vec2 e = v[j] - v[i];
        vec2 w = p - v[i];
        vec2 b = w - e * clamp(dot(w, e) / dot(e, e), 0.0f, 1.0f);
        d = min(d, dot(b, b));
        bvec3 c = bvec3(p.y >= v[i].y, p.y<v[j].y, e.x* w.y>e.y * w.x);
        if (all(c) || all(not(c)))
            s *= -1.0;
    }
    return s * sqrt(d);
}

bind(sdEllipse(std::placeholders::_1, vec2 ab)
{
    p = abs(p);
    if (p.x > p.y)
    {
        p = vec2(p.y, p.x);
        ab = vec2(ab.y, ab.x);
    }
    float l = ab.y * ab.y - ab.x * ab.x;
    float m = ab.x * p.x / l;      float m2 = m * m;
    float n = ab.y * p.y / l;      float n2 = n * n;
    float c = (m2 + n2 - 1.0) / 3.0; float c3 = c * c * c;
    float q = c3 + m2 * n2 * 2.0;
    float d = c3 + m2 * n2;
    float g = m + m * n2;
    float co;
    if (d < 0.0)
    {
        float h = acos(q / c3) / 3.0;
        float s = cos(h);
        float t = sin(h) * sqrt(3.0);
        float rx = sqrt(-c * (s + t + 2.0) + m2);
        float ry = sqrt(-c * (s - t + 2.0) + m2);
        co = (ry + sign(l) * rx + abs(g) / (rx * ry) - m) / 2.0;
    }
    else
    {
        float h = 2.0 * m * n * sqrt(d);
        float s = sign(q + h) * pow(abs(q + h), 1.0 / 3.0);
        float u = sign(q - h) * pow(abs(q - h), 1.0 / 3.0);
        float rx = -s - u - c * 4.0 + 2.0 * m2;
        float ry = (s - u) * sqrt(3.0);
        float rm = sqrt(rx * rx + ry * ry);
        co = (ry / sqrt(rm - rx) + 2.0 * g / rm - m) / 2.0;
    }
    vec2 r = ab * vec2(co, sqrt(1.0 - co * co));
    return length(r - p) * sign(p.y - r.y);
}

bind(sdParabola(vec2 pos, float k)
{
    pos.x = abs(pos.x);
    float ik = 1.0 / k;
    float p = ik * (pos.y - 0.5 * ik) / 3.0;
    float q = 0.25 * ik * ik * pos.x;
    float h = q * q - p * p * p;
    float r = sqrt(abs(h));
    float x = (h > 0.0) ?
        pow(q + r, 1.0 / 3.0) - pow(abs(q - r), 1.0 / 3.0) * sign(r - q) :
        2.0 * cos(atan2(r, q) / 3.0) * sqrt(p);
    return length(pos - vec2(x, k * x * x)) * sign(pos.x - x);
}

bind(sdParabola(vec2 pos, float wi, float he)
{
    pos.x = abs(pos.x);
    float ik = wi * wi / he;
    float p = ik * (he - pos.y - 0.5 * ik) / 3.0;
    float q = pos.x * ik * ik * 0.25;
    float h = q * q - p * p * p;
    float r = sqrt(abs(h));
    float x = (h > 0.0) ?
        pow(q + r, 1.0 / 3.0) - pow(abs(q - r), 1.0 / 3.0) * sign(r - q) :
        2.0 * cos(atan(r / q) / 3.0) * sqrt(p);
    x = min(x, wi);
    return length(pos - vec2(x, he - x * x / ik)) *
        sign(ik * (pos.y - he) + pos.x * pos.x);
}

bind(sdBezier(vec2 pos, vec2 A, vec2 B, vec2 C)
{
    vec2 a = B - A;
    vec2 b = A - 2.0 * B + C;
    vec2 c = a * 2.0;
    vec2 d = A - pos;
    float kk = 1.0 / dot(b, b);
    float kx = kk * dot(a, b);
    float ky = kk * (2.0 * dot(a, a) + dot(d, b)) / 3.0;
    float kz = kk * dot(d, a);
    float res = 0.0;
    float p = ky - kx * kx;
    float p3 = p * p * p;
    float q = kx * (2.0 * kx * kx - 3.0 * ky) + kz;
    float h = q * q + 4.0 * p3;
    if (h >= 0.0)
    {
        h = sqrt(h);
        vec2 x = (vec2(h, -h) - q) / 2.0;
        vec2 uv = sign(x) * pow(abs(x), vec2(1.0f / 3.0f));
        float t = clamp(uv.x + uv.y - kx, 0.0f, 1.0f);
        res = dot2(d + (c + b * t) * t);
    }
    else
    {
        float z = sqrt(-p);
        float v = acos(q / (p * z * 2.0)) / 3.0;
        float m = cos(v);
        float n = sin(v) * 1.732050808;
        vec3  t = clamp(vec3(m + m, -n - m, n - m) * z - kx, 0.0, 1.0);
        res = min(dot2(d + (c + b * t.x) * t.x),
            dot2(d + (c + b * t.y) * t.y));
        // the third root cannot be the closest
        // res = min(res,dot2(d+(c+b*t.z)*t.z));
    }
    return sqrt(res);
}
*/