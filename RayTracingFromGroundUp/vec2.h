#ifndef _vec2_h_
#define _vec2_h_

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

#endif