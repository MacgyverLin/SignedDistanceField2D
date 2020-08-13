#ifndef _vec4_h_
#define _vec4_h_

#include "util.h"

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

#endif