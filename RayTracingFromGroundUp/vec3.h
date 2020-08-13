#ifndef _vec3_h_
#define _vec3_h_

#include "util.h"

//////////////////////////////////////////////////////////
class vec3
{
public:
	vec3(float v = 0)
	{
		x = v;
		y = v;
		z = v;
	}

	vec3(float e0, float e1, float e2)
	{
		x = e0;
		y = e1;
		z = e2;
	}

	vec3(const vec3& other)
	{
		x = other.x;
		y = other.y;
		z = other.z;
	}

	vec3& operator = (const vec3& other)
	{
		x = other.x;
		y = other.y;
		z = other.z;

		return *this;
	}

	inline float operator[] (int i) const { return ((float*)(&x))[i]; }
	inline float& operator[] (int i) { return ((float*)(&x))[i]; }

	inline vec3 operator+() const { return *this; }
	inline vec3 operator-() const { return vec3(-x, -y, -z); }

	inline vec3& operator += (const vec3& other)
	{
		x += other.x;
		y += other.y;
		z += other.z;
		return *this;
	}

	inline vec3& operator -= (const vec3& other)
	{
		x -= other.x;
		y -= other.y;
		z -= other.z;
		return *this;
	}

	inline vec3& operator *= (const vec3& other)
	{
		x *= other.x;
		y *= other.y;
		z *= other.z;
		return *this;
	}

	inline vec3& operator /= (const vec3& other)
	{
		x /= other.x;
		y /= other.y;
		z /= other.z;
		return *this;
	}

	inline vec3& operator *= (const float t)
	{
		x *= t;
		y *= t;
		z *= t;
		return *this;
	}

	inline vec3& operator /= (const float t)
	{
		x /= t;
		y /= t;
		z /= t;
		return *this;
	}

	friend inline vec3 operator + (const vec3& v1, const vec3& v2)
	{
		return vec3(v1.x + v2.x, v1.y + v2.y, v1.z + v2.z);
	}

	friend inline vec3 operator - (const vec3& v1, const vec3& v2)
	{
		return vec3(v1.x - v2.x, v1.y - v2.y, v1.z - v2.z);
	}

	friend inline vec3 operator * (const vec3& v1, const vec3& v2)
	{
		return vec3(v1.x * v2.x, v1.y * v2.y, v1.z * v2.z);
	}

	friend inline vec3 operator / (const vec3& v1, const vec3& v2)
	{
		return vec3(v1.x / v2.x, v1.y / v2.y, v1.z / v2.z);
	}

	friend inline vec3 operator * (const vec3& v1, const float& t)
	{
		return vec3(v1.x * t, v1.y * t, v1.z * t);
	}

	friend inline vec3 operator * (const float& t, const vec3& v1)
	{
		return vec3(v1.x * t, v1.y * t, v1.z * t);
	}

	friend inline vec3 operator / (const vec3& v1, const float& t)
	{
		return vec3(v1.x / t, v1.y / t, v1.z / t);
	}

	inline float length() const
	{
		return sqrt(square_length());
	}

	inline float square_length() const
	{
		return x * x + y * y + z * z;
	}

	inline void make_unit_vector()
	{
		*this = *this / length();
	}

	friend inline float dot(const vec3& v1, const vec3& v2)
	{
		return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
	}

	friend inline vec3 cross(const vec3& v1, const vec3& v2)
	{
		return vec3
		(
			(v1.y * v2.z - v1.z * v2.y),
			-(v1.x * v2.z - v1.z * v2.x),
			(v1.x * v2.y - v1.y * v2.x)
		);
	}

	friend inline vec3 unit_vector(vec3 v)
	{
		float l = v.length();
		return v / l;
	}

	friend inline vec3 abs(const vec3& v)
	{
		return vec3
		(
			abs(v.x),
			abs(v.y),
			abs(v.z)
		);
	}

	friend inline vec3 max(const vec3& v1, const vec3& v2)
	{
		return vec3
		(
			max(v1.x, v2.x),
			max(v1.y, v2.y),
			max(v1.z, v2.z)
		);
	}

	friend inline vec3 min(const vec3& v1, const vec3& v2)
	{
		return vec3
		(
			min(v1.x, v2.x),
			min(v1.y, v2.y),
			max(v1.z, v2.z)
		);
	}

	friend inline vec3 mix(const vec3& v1, const vec3& v2, float t)
	{
		return v1 * (1 - t) + v2 * t;
	}

	friend inline vec3 sin(const vec3& v)
	{
		return vec3
		(
			sin(v.x),
			sin(v.y),
			sin(v.z)
		);
	}

	friend inline vec3 cos(const vec3& v)
	{
		return vec3
		(
			cos(v.x),
			cos(v.y),
			cos(v.z)
		);
	}

	friend inline vec3 tan(const vec3& v)
	{
		return vec3
		(
			tan(v.x),
			tan(v.y),
			tan(v.z)
		);
	}

	friend inline vec3 mod(const vec3& v, const vec3& c)
	{
		return vec3
		(
			fmod(v.x, c.x),
			fmod(v.y, c.y),
			fmod(v.z, c.z)
		);
	}

	float x;
	float y;
	float z;
};

#endif