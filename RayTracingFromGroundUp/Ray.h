#ifndef _RAY_h_
#define _RAY_h_

#include "util.h"

class vec2
{
public:
	vec2(float v = 0)
	{
		x = v;
		y = v;
	}

	vec2(float e0, float e1)
	{
		x = e0;
		y = e1;
	}

	vec2(const vec2& other)
	{
		x = other.x;
		y = other.y;
	}

	vec2& operator = (const vec2& other)
	{
		x = other.x;
		y = other.y;

		return *this;
	}

	inline float operator[] (int i) const { return ((float*)(&x))[i]; }
	inline float& operator[] (int i) { return ((float*)(&x))[i]; }

	inline vec2 operator+() const { return *this; }
	inline vec2 operator-() const { return vec2(-x, -y); }

	inline vec2& operator += (const vec2& other)
	{
		x += other.x;
		y += other.y;
		return *this;
	}

	inline vec2& operator -= (const vec2& other)
	{
		x -= other.x;
		y -= other.y;
		return *this;
	}

	inline vec2& operator *= (const vec2& other)
	{
		x *= other.x;
		y *= other.y;
		return *this;
	}

	inline vec2& operator /= (const vec2& other)
	{
		x /= other.x;
		y /= other.y;
		return *this;
	}

	inline vec2& operator *= (const float t)
	{
		x *= t;
		y *= t;
		return *this;
	}

	inline vec2& operator /= (const float t)
	{
		x /= t;
		y /= t;
		return *this;
	}

	friend inline vec2 operator + (const vec2& v1, const vec2& v2)
	{
		return vec2(v1.x + v2.x, v1.y + v2.y);
	}

	friend inline vec2 operator - (const vec2& v1, const vec2& v2)
	{
		return vec2(v1.x - v2.x, v1.y - v2.y);
	}

	friend inline vec2 operator * (const vec2& v1, const vec2& v2)
	{
		return vec2(v1.x * v2.x, v1.y * v2.y);
	}

	friend inline vec2 operator / (const vec2& v1, const vec2& v2)
	{
		return vec2(v1.x / v2.x, v1.y / v2.y);
	}

	friend inline vec2 operator * (const vec2& v1, const float& t)
	{
		return vec2(v1.x * t, v1.y * t);
	}

	friend inline vec2 operator * (const float& t, const vec2& v1)
	{
		return vec2(v1.x * t, v1.y * t);
	}

	friend inline vec2 operator / (const vec2& v1, const float& t)
	{
		return vec2(v1.x / t, v1.y / t);
	}

	inline void make_zero()
	{
		x = 0;
		y = 0;
	}

	inline void make_unitx()
	{
		x = 1;
		y = 0;
	}

	inline void make_unity()
	{
		x = 0;
		y = 1;
	}

	inline float length() const 
	{ 
		return sqrt(square_length()); 
	}

	inline float square_length() const 
	{ 
		return x * x + y * y;
	}

	inline vec2 perp() const
	{
		return vec2(-x, y);
	}

	inline const vec2& normalize()
	{
		*this = *this / length();

		return *this;
	}

	inline float dot(const vec2& v) const
	{
		return x * v.x + y * v.y;
	}

	inline float ndot(const vec2& v) const
	{
		return x * v.x - y * v.y;
	}

	friend inline float square_length(const vec2& v)
	{
		return v.square_length();
	}

	friend inline float length(const vec2& v)
	{
		return v.length();
	}

	friend inline float dot(const vec2& v1, const vec2& v2)
	{
		return v1.dot(v2);
	}

	friend inline float ndot(const vec2& v1, const vec2& v2)
	{
		return v1.ndot(v2);
	}

	friend inline vec2 pow(const vec2& v, const vec2& p)
	{
		return vec2(pow(v.x, p.x), pow(v.y, p.y));
	}

	friend inline vec2 sign(const vec2& v)
	{
		return vec2(sign(v.x), sign(v.y));
	}

	friend inline vec2 abs(const vec2& v)
	{
		return vec2(abs(v.x), abs(v.y));
	}

	friend inline vec2 max(const vec2& v1, const vec2& v2)
	{
		return vec2
		(
			max(v1.x, v2.x), 
			max(v1.y, v2.y)
		);
	}

	friend inline vec2 min(const vec2& v1, const vec2& v2)
	{
		return vec2
		(
			min(v1.x, v2.x), 
			min(v1.y, v2.y)
		);
	}

	friend inline vec2 mix(const vec2& v1, const vec2& v2, float t)
	{
		return v1 * (1 - t) + v2;
	}

	friend inline vec2 sin(const vec2& v)
	{
		return vec2
		(
			sin(v.x),
			sin(v.y)
		);
	}

	friend inline vec2 cos(const vec2& v)
	{
		return vec2
		(
			cos(v.x),
			cos(v.y)
		);
	}

	friend inline vec2 tan(const vec2& v)
	{
		return vec2
		(
			tan(v.x),
			tan(v.y)
		);
	}

	friend inline vec2 mod(const vec2& v, const vec2& c)
	{
		return vec2
		(
			fmod(v.x, c.x),
			fmod(v.y, c.y)
		);
	}

	float x;
	float y;
};

//////////////////////////////////////////////////////////
class vec4
{
public:
	vec4(float v = 0)
	{
		x = v;
		y = v;
		z = v;
		w = v;
	}

	vec4(float e0, float e1, float e2, float e3)
	{
		x = e0;
		y = e1;
		z = e2;
		w = e3;
	}

	vec4(const vec4& other)
	{
		x = other.x;
		y = other.y;
		z = other.z;
		w = other.w;
	}

	vec4& operator = (const vec4& other)
	{
		x = other.x;
		y = other.y;
		z = other.z;
		w = other.w;

		return *this;
	}

	inline float operator[] (int i) const { return ((float*)(&x))[i]; }
	inline float& operator[] (int i) { return ((float*)(&x))[i]; }

	inline vec4 operator+() const { return *this; }
	inline vec4 operator-() const { return vec4(-x, -y, -z, -w); }

	inline vec4& operator += (const vec4& other)
	{
		x += other.x;
		y += other.y;
		z += other.z;
		w += other.w;
		return *this;
	}

	inline vec4& operator -= (const vec4& other)
	{
		x -= other.x;
		y -= other.y;
		z -= other.z;
		w -= other.w;
		return *this;
	}

	inline vec4& operator *= (const vec4& other)
	{
		x *= other.x;
		y *= other.y;
		z *= other.z;
		w *= other.w;
		return *this;
	}

	inline vec4& operator /= (const vec4& other)
	{
		x /= other.x;
		y /= other.y;
		z /= other.z;
		w /= other.w;
		return *this;
	}

	inline vec4& operator *= (const float t)
	{
		x *= t;
		y *= t;
		z *= t;
		w *= t;
		return *this;
	}

	inline vec4& operator /= (const float t)
	{
		x /= t;
		y /= t;
		z /= t;
		w /= t;
		return *this;
	}

	friend inline vec4 operator + (const vec4& v1, const vec4& v2)
	{
		return vec4(v1.x + v2.x, v1.y + v2.y, v1.z + v2.z, v1.w + v2.w);
	}

	friend inline vec4 operator - (const vec4& v1, const vec4& v2)
	{
		return vec4(v1.x - v2.x, v1.y - v2.y, v1.z - v2.z, v1.w - v2.w);
	}

	friend inline vec4 operator * (const vec4& v1, const vec4& v2)
	{
		return vec4(v1.x * v2.x, v1.y * v2.y, v1.z * v2.z, v1.w * v2.w);
	}

	friend inline vec4 operator / (const vec4& v1, const vec4& v2)
	{
		return vec4(v1.x / v2.x, v1.y / v2.y, v1.z / v2.z, v1.w / v2.w);
	}

	friend inline vec4 operator * (const vec4& v1, const float& t)
	{
		return vec4(v1.x * t, v1.y * t, v1.z * t, v1.w * t);
	}

	friend inline vec4 operator * (const float& t, const vec4& v1)
	{
		return vec4(v1.x * t, v1.y * t, v1.z * t, v1.w * t);
	}

	friend inline vec4 operator / (const vec4& v1, const float& t)
	{
		return vec4(v1.x / t, v1.y / t, v1.z / t, v1.w / t);
	}

	inline void make_zero()
	{
		x = 0;
		y = 0;
		z = 0;
		w = 0;
	}

	inline void make_unitx()
	{
		x = 1;
		y = 0;
		z = 0;
		w = 0;
	}

	inline void make_unity()
	{
		x = 0;
		y = 1;
		z = 0;
		w = 0;
	}

	inline void make_unitz()
	{
		x = 0;
		y = 0;
		z = 1;
		w = 0;
	}

	inline void make_unitw()
	{
		x = 0;
		y = 0;
		z = 0;
		w = 1;
	}

	inline float length() const
	{
		return sqrt(square_length());
	}

	inline float square_length() const
	{
		return x * x + y * y + z * z + w * w;
	}

	inline const vec4& normalize()
	{
		*this = *this / length();

		return *this;
	}

	inline float dot(const vec4& v) const
	{
		return x * v.x + y * v.y + z * v.z + w * v.w;
	}

	friend inline float square_length(const vec4& v)
	{
		return v.square_length();
	}

	friend inline float length(const vec4& v)
	{
		return v.length();
	}

	friend inline float dot(const vec4& v1, const vec4& v2)
	{
		return v1.dot(v2);
	}

	friend inline vec4 abs(const vec4& v)
	{
		return vec4
		(
			abs(v.x), 
			abs(v.y), 
			abs(v.z), 
			abs(v.w)
		);
	}

	friend inline vec4 max(const vec4& v1, const vec4& v2)
	{
		return vec4
		(
			max(v1.x, v2.x), 
			max(v1.y, v2.y), 
			max(v1.z, v2.z), 
			max(v1.w, v2.w)
		);
	}

	friend inline vec4 min(const vec4& v1, const vec4& v2)
	{
		return vec4
		(
			min(v1.x, v2.x), 
			min(v1.y, v2.y), 
			min(v1.z, v2.z), 
			min(v1.w, v2.w)
		);
	}

	friend inline vec4 mix(const vec4& v1, const vec4& v2, float t)
	{
		return v1 * (1 - t) + v2;
	}

	friend inline vec4 sin(const vec4& v)
	{
		return vec4
		(
			sin(v.x),
			sin(v.y),
			sin(v.z),
			sin(v.w)
		);
	}

	friend inline vec4 cos(const vec4& v)
	{
		return vec4
		(
			cos(v.x),
			cos(v.y),
			cos(v.z),
			cos(v.w)
		);
	}

	friend inline vec4 tan(const vec4& v)
	{
		return vec4
		(
			tan(v.x),
			tan(v.y),
			tan(v.z),
			tan(v.w)
		);
	}

	friend inline vec4 mod(const vec4& v, const vec4& c)
	{
		return vec4
		(
			fmod(v.x, c.x),
			fmod(v.y, c.y),
			fmod(v.z, c.z),
			fmod(v.w, c.w)
		);
	}

	float x;
	float y;
	float z;
	float w;
};

//////////////////////////////////////////////////////////
class mat2
{
public:
	mat2(float e00, float e01, float e10, float e11)
	{
		e[0][0] = e00;
		e[0][1] = e01;
		e[1][0] = e10;
		e[1][1] = e11;
	}

	mat2(const mat2& other)
	{
		e[0][0] = other.e[0][0];
		e[0][1] = other.e[0][1];
		e[1][0] = other.e[1][0];
		e[1][1] = other.e[1][1];
	}

	mat2& operator = (const mat2& other)
	{
		e[0][0] = other.e[0][0];
		e[0][1] = other.e[0][1];
		e[1][0] = other.e[1][0];
		e[1][1] = other.e[1][1];

		return *this;
	}

	inline mat2 operator+() const { return *this; }
	inline mat2 operator-() const { return mat2(-e[0][0], -e[0][1], -e[1][0], -e[1][1]); }

	inline mat2& operator += (const mat2& other)
	{
		e[0][0] += other.e[0][0];
		e[0][1] += other.e[0][1];
		e[1][0] += other.e[1][0];
		e[1][1] += other.e[1][1];

		return *this;
	}

	inline mat2& operator -= (const mat2& other)
	{
		e[0][0] -= other.e[0][0];
		e[0][1] -= other.e[0][1];
		e[1][0] -= other.e[1][0];
		e[1][1] -= other.e[1][1];

		return *this;
	}

	inline mat2& operator *= (const mat2& other)
	{
		e[0][0] *= other.e[0][0];
		e[0][1] *= other.e[0][1];
		e[1][0] *= other.e[1][0];
		e[1][1] *= other.e[1][1];

		return *this;
	}

	inline mat2& operator /= (const mat2& other)
	{
		e[0][0] /= other.e[0][0];
		e[0][1] /= other.e[0][1];
		e[1][0] /= other.e[1][0];
		e[1][1] /= other.e[1][1];

		return *this;
	}

	inline mat2& operator *= (const float t)
	{
		e[0][0] *= t;
		e[0][1] *= t;
		e[1][0] *= t;
		e[1][1] *= t;

		return *this;
	}

	inline mat2& operator /= (const float t)
	{
		e[0][0] /= t;
		e[0][1] /= t;
		e[1][0] /= t;
		e[1][1] /= t;

		return *this;
	}

	friend inline mat2 operator + (const mat2& m0, const mat2& m1)
	{
		return mat2(m0.e[0][0] + m1.e[0][0], m0.e[0][1] + m1.e[0][1], m0.e[1][0] + m1.e[1][0], m0.e[1][1] + m1.e[1][1]);
	}

	friend inline mat2 operator - (const mat2& m0, const mat2& m1)
	{
		return mat2(m0.e[0][0] - m1.e[0][0], m0.e[0][1] - m1.e[0][1], m0.e[1][0] - m1.e[1][0], m0.e[1][1] - m1.e[1][1]);
	}

	friend inline mat2 operator * (const mat2& m0, const mat2& m1)
	{
		mat2 r(0, 0, 0, 0);
		
		r.e[0][0] = m0.e[0][0] * m1.e[0][0] + m0.e[0][1] * m1.e[1][0];
		r.e[0][1] = m0.e[0][0] * m1.e[0][1] + m0.e[0][1] * m1.e[1][1];
		r.e[1][0] = m0.e[1][0] * m1.e[0][0] + m0.e[1][1] * m1.e[1][0];
		r.e[1][1] = m0.e[1][0] * m1.e[0][1] + m0.e[1][1] * m1.e[1][1];

		return r;
	}

	friend inline vec2 operator * (const mat2& m0, const vec2& v)
	{
		vec2 r(0, 0);

		r.x = m0.e[0][0] * v.x + m0.e[0][1] * v.y;
		r.y = m0.e[0][0] * v.x + m0.e[0][1] * v.y;

		return r;
	}

	friend inline mat2 operator * (const mat2& m0, const float& t)
	{
		return mat2(m0.e[0][0] * t, m0.e[0][1] * t, m0.e[1][0] * t, m0.e[1][1] * t);
	}

	friend inline mat2 operator * (const float& t, const mat2& m0)
	{
		return mat2(m0.e[0][0] * t, m0.e[0][1] * t, m0.e[1][0] * t, m0.e[1][1] * t);
	}

	friend inline mat2 operator / (const mat2& m0, const float& t)
	{
		return mat2(m0.e[0][0] / t, m0.e[0][1] / t, m0.e[1][0] / t, m0.e[1][1] / t);
	}

	inline void make_zero()
	{
		e[0][0] = 0;
		e[0][1] = 0;
		e[1][0] = 0;
		e[1][1] = 0;
	}

	inline void make_identity()
	{
		e[0][0] = 1;
		e[0][1] = 0;
		e[1][0] = 0;
		e[1][1] = 1;
	}

	inline void make_rotate(float angle)
	{
		float s = sin(angle / 180 * 3.1416);
		float c = cos(angle / 180 * 3.1416);

		e[0][0] = c;
		e[0][1] = -s;
		e[1][0] = s;
		e[1][1] = c;
	}

	inline float determinant() const
	{
		return (e[0][0] * e[1][1]) - (e[0][1] * e[1][0]);
	}

	inline mat2 inverse() const
	{
		float d = determinant();
		return mat2(e[1][1], -e[0][1], -e[1][0], e[0][0]) /d;
	}

	friend inline mat2 inverse(const mat2& m)
	{
		return m.inverse();
	}
	
	friend inline mat2 abs(const mat2& m0)
	{
		return mat2
		(
			abs(m0.e[0][0]), 
			abs(m0.e[0][1]), 
			abs(m0.e[1][0]), 
			abs(m0.e[1][1])
		);
	}

	friend inline mat2 max(const mat2& m0, const mat2& m1)
	{
		return mat2
		(
			max(m0.e[0][0], m1.e[0][0]),
			max(m0.e[0][1], m1.e[0][1]),
			max(m0.e[1][0], m1.e[1][0]),
			max(m0.e[1][1], m1.e[1][1])
		);
	}

	friend inline mat2 min(const mat2& m0, const mat2& m1)
	{
		return mat2
		(
			min(m0.e[0][0], m1.e[0][0]),
			min(m0.e[0][1], m1.e[0][1]),
			min(m0.e[1][0], m1.e[1][0]),
			min(m0.e[1][1], m1.e[1][1])
		);
	}

	friend inline mat2 mix(const mat2& m0, const mat2& m1, float t)
	{
		return m0 * (1 - t) + m1;
	}

	float e[2][2];
};

#endif