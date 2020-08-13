#ifndef _mat2_h_
#define _mat2_h_

#include "util.h"
#include "vec2.h"

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