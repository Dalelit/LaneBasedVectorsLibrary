#pragma once

///////////////////////////////////////////////////
// TO DO
//
// - work out constants... for zero and one

#define _USE_MATH_DEFINES 1

// C RunTime Header Files
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <malloc.h>
#include <memory.h>
#include <tchar.h>
#include <math.h>
#include <time.h>
#include <emmintrin.h>

#define LANE_WIDTH 4
#define LANE_WIDTH_FLOAT 4.0f;

#if LANE_WIDTH == 1

struct LaneFloat
{
	float value;

	LaneFloat(float _value = 0.0f) { value = _value; }

	LaneFloat operator+  (const LaneFloat &b) { return value + b.value; }
	LaneFloat operator-  (const LaneFloat &b) { return value - b.value; }
	LaneFloat operator*  (const LaneFloat &b) { return value * b.value; }
	LaneFloat operator/  (const LaneFloat &b) { return value / b.value; }
	LaneFloat operator+= (const LaneFloat &b) { return value += b.value; }
	LaneFloat operator-= (const LaneFloat &b) { return value -= b.value; }
	LaneFloat operator*= (const LaneFloat &b) { return value *= b.value; }

	LaneFloat operator-  () { return -value; }

	bool operator>  (const LaneFloat &b) { return value > b.value; }
	bool operator>= (const LaneFloat &b) { return value >= b.value; }
	bool operator<  (const LaneFloat &b) { return value < b.value; }
	bool operator<= (const LaneFloat &b) { return value <= b.value; }
	bool operator>  (float b) { return value > b; }
	bool operator>= (float b) { return value >= b; }
	bool operator<  (float b) { return value < b; }
	bool operator<= (float b) { return value <= b; }

	LaneFloat sqrt() { return sqrtf(value); }
	LaneFloat abs() { return LaneFloat(fabsf(value)); }
};


#elif LANE_WIDTH == 4


struct LaneFloat
{
	__m128 value;

	LaneFloat() {};
	LaneFloat(__m128 _value) { value = _value; }
	LaneFloat(float f) { *this = f; };

	void Load(float* f) { value = _mm_load_ps(f); }
	void Store(float* f) { _mm_store_ps(f, value); }

	LaneFloat operator= (const float &b) { value = _mm_load_ps1(&b); return *this; }
	LaneFloat operator= (const LaneFloat &b) { value = b.value; return *this; }

	LaneFloat operator+   (const LaneFloat &b) { return _mm_add_ps(value, b.value); }
	LaneFloat operator-   (const LaneFloat &b) { return _mm_sub_ps(value, b.value); }
	LaneFloat operator*   (const LaneFloat &b) { return _mm_mul_ps(value, b.value); }
	LaneFloat operator/   (const LaneFloat &b) { return _mm_div_ps(value, b.value); }
	LaneFloat operator+=  (const LaneFloat &b) { value = _mm_add_ps(value, b.value); return *this; }
	LaneFloat operator-=  (const LaneFloat &b) { value = _mm_sub_ps(value, b.value); return *this; }
	LaneFloat operator*=  (const LaneFloat &b) { value = _mm_mul_ps(value, b.value); return *this; }
	LaneFloat operator/=  (const LaneFloat &b) { value = _mm_div_ps(value, b.value); return *this; }

	LaneFloat operator-  () { return _mm_sub_ps(_mm_setzero_ps(), value); }

	LaneFloat operator== (const LaneFloat &b) { return _mm_cmpeq_ps(value, b.value); }
	LaneFloat operator!= (const LaneFloat &b) { return _mm_cmpneq_ps(value, b.value); }
	LaneFloat operator>  (const LaneFloat &b) { return _mm_cmpgt_ps(value, b.value); }
	LaneFloat operator>= (const LaneFloat &b) { return _mm_cmpge_ps(value, b.value); }
	LaneFloat operator<  (const LaneFloat &b) { return _mm_cmplt_ps(value, b.value); }
	LaneFloat operator<= (const LaneFloat &b) { return _mm_cmple_ps(value, b.value); }

	float HorizAdd() { return value.m128_f32[0] + value.m128_f32[1] + value.m128_f32[2] + value.m128_f32[3]; } // TO DO - is there a better way?
	float HorizAvg() { return HorizAdd() / LANE_WIDTH_FLOAT; }

	bool MaskIsZero() { return HorizAdd() == 0.0f; }

	LaneFloat MaskValue(const LaneFloat &mask) { return _mm_and_ps(mask.value, value); }

	LaneFloat MaskedAssign(const LaneFloat &mask, const LaneFloat &b) { value = _mm_add_ps(_mm_andnot_ps(mask.value, value), _mm_and_ps(mask.value, b.value)); return *this; }

	LaneFloat sqrt() { return _mm_sqrt_ps(value); }

	LaneFloat abs() { // TO DO: this feels like a hack
		__m128 mask = _mm_cmpge_ps(value, _mm_setzero_ps());
		__m128 negated = _mm_andnot_ps(mask, value); // only the ones that are negative
		return _mm_sub_ps(_mm_sub_ps(value, negated), negated); // add them twice to make positive
		;
	}

	LaneFloat clamp(LaneFloat min, LaneFloat max) { MaskedAssign(*this < min, min); MaskedAssign(*this > max, max); return *this; }
};

#endif

struct Vector
{
	float x;
	float y;
	float z;

	Vector(float _x = 0, float _y = 0, float _z = 0) { x = _x; y = _y; z = _z; }
};

struct LaneVector
{
	LaneFloat x;
	LaneFloat y;
	LaneFloat z;

	LaneVector(LaneFloat _x = 0, LaneFloat _y = 0, LaneFloat _z = 0) { x = _x; y = _y; z = _z; }

	void Load(Vector* v) { // TO DO: should this use the LaneFloat store function?
		x.value.m128_f32[0] = v[0].x; y.value.m128_f32[0] = v[0].y; z.value.m128_f32[0] = v[0].z;
		x.value.m128_f32[1] = v[1].x; y.value.m128_f32[1] = v[1].y; z.value.m128_f32[1] = v[1].z;
		x.value.m128_f32[2] = v[2].x; y.value.m128_f32[2] = v[2].y; z.value.m128_f32[2] = v[2].z;
		x.value.m128_f32[3] = v[3].x; y.value.m128_f32[3] = v[3].y; z.value.m128_f32[3] = v[3].z;
	}

	void Store(Vector* v) {
		v[0].x = x.value.m128_f32[0]; v[0].y = y.value.m128_f32[0]; v[0].z = z.value.m128_f32[0];
		v[1].x = x.value.m128_f32[1]; v[1].y = y.value.m128_f32[1]; v[1].z = z.value.m128_f32[1];
		v[2].x = x.value.m128_f32[2]; v[2].y = y.value.m128_f32[2]; v[2].z = z.value.m128_f32[2];
		v[3].x = x.value.m128_f32[3]; v[3].y = y.value.m128_f32[3]; v[3].z = z.value.m128_f32[3];
	}

	// TO DO: have the operators with the floats on the left

	LaneVector operator=(const LaneFloat& b) { x = b; y = b; z = b; return *this; }
	//LaneVector operator= (const LaneVector &b) { x = b.x; y = b.y, z = b.z; return *this; }
	LaneVector operator+ (const LaneVector &b) { return LaneVector(x + b.x, y + b.y, z + b.z); }
	LaneVector operator- (const LaneVector &b) { return LaneVector(x - b.x, y - b.y, z - b.z); }
	LaneVector operator* (const LaneVector &b) { return LaneVector(x * b.x, y * b.y, z * b.z); }
	LaneVector operator* (const LaneFloat &b) { return LaneVector(x * b, y * b, z * b); }
	LaneVector operator/ (const LaneFloat &b) { return LaneVector(x / b, y / b, z / b); }
	LaneVector operator+= (const LaneVector &b) { x += b.x; y += b.y, z += b.z; return *this; }
	LaneVector operator*= (const LaneFloat &b) { x *= b; y *= b, z *= b; return *this; }

	LaneFloat  dot(const LaneVector &b) { return x*b.x + y*b.y + z*b.z; }
	LaneVector cross(const LaneVector &b) { return LaneVector(y*b.z - z*b.y, z*b.x - x*b.z, x*b.y - y*b.x); }

	LaneFloat  lengthSqrd() { return x*x + y*y + z*z; }
	LaneFloat  length() { return (x*x + y*y + z*z).sqrt(); }
	LaneVector normalize() { return *this = *this / this->length(); }

	Vector HorizAdd() { return Vector(x.HorizAdd(), y.HorizAdd(), z.HorizAdd()); }
	Vector HorizAvg() { return Vector(x.HorizAvg(), y.HorizAvg(), z.HorizAvg()); }

	LaneVector clamp(LaneFloat min, LaneFloat max) { x.clamp(min, max); y.clamp(min, max); z.clamp(min, max);
	}
};

struct LaneUInt
{
	__m128i value;

	LaneUInt() {}
	LaneUInt(__m128i _value) { value = _value; }
	LaneUInt(int i) { value.m128i_i32[0] = i; value.m128i_i32[1] = i; value.m128i_i32[2] = i; value.m128i_i32[3] = i; }

	LaneUInt operator<< (const int &b) { return _mm_slli_epi32(value, b); }
	LaneUInt operator>> (const int &b) { return _mm_srli_epi32(value, b); }

	LaneUInt operator^  (const LaneUInt &b) { return _mm_xor_si128(value, b.value); }
	LaneUInt operator^= (const LaneUInt &b) { value = _mm_xor_si128(value, b.value); return *this; }

	LaneFloat ToFloat() { return LaneFloat(_mm_cvtepi32_ps(value)); }
};

LaneUInt LaneSeed();
float xorshift32(UINT32* state);
LaneFloat LaneXorShift32(LaneUInt* state);




// TO DO: sort these out properly
inline float RandomUnilateral()
{
	return (float)rand() / (float)RAND_MAX;
}

// TO DO: sort these out properly
inline float RandomBilateral()
{
	return -1.0f + 2.0f * RandomUnilateral();
}

// TO DO: sort these out properly
inline UINT32 Random0To255()
{
	return (UINT32)(255.0f * RandomUnilateral());
}
