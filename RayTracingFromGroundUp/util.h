#ifndef _System_h_
#define _System_h_

#include <string.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <float.h>
#include <math.h>
#include <memory>
#include <vector>
#include <map>
#include <algorithm>
#include <functional>
using namespace std;

//#define ffmax(a, b) (((a)>(b)) ? (a) : (b))
//#define ffmin(a, b) (((a)<(b)) ? (a) : (b))
#define clamp(x, a, b) (max(min((x), (b)), (a)))
//#define sign(a) (((a)>0)?1:(((a)<0)?-1:0))

float mix(float a, float b, float t)
{
	return a * (1 - t) + b * t;
}

template <typename T> int sign(T val) {
	return (T(0) < val) - (val < T(0));
}

template<class T>
float random(const T& min, const T& max)
{
	return rand() * (max - min) / RAND_MAX + min;
}

template<class T>
T smoothstep(const T& edge0, const T& edge1, const T& x)
{
	T t = clamp((x - edge0) / (edge1 - edge0), 0.0f, 1.0f);
	return (3.0 - 2.0 * t) * t * t;
}

#endif