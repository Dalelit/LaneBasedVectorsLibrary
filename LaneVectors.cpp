#include "stdafx.h"
#include "LaneVectors.h"

LaneUInt LaneSeed()
{
	LaneUInt seed;

	seed.value.m128i_u32[0] = rand();
	seed.value.m128i_u32[1] = rand();
	seed.value.m128i_u32[2] = rand();
	seed.value.m128i_u32[3] = rand();

	return seed;
}

float xorshift32(UINT32* state)
{
	// Algorithm "xor" from p. 4 of Marsaglia, "Xorshift RNGs"
	UINT32 x = state[0];
	x ^= x << 13;
	x ^= x >> 17;
	x ^= x << 5;
	state[0] = x;
	return (float)x / (float)UINT32_MAX;
}

LaneFloat LaneXorShift32(LaneUInt* state)
{
	LaneUInt x = *state;
	x ^= x << 13;
	x ^= x >> 17;
	x ^= x << 5;
	*state = x;

	LaneFloat a = x.ToFloat();
	LaneUInt ui = LaneUInt(INT32_MAX - 1);
	LaneFloat b = ui.ToFloat();

	return a / b;
}


/*
inline void printLaneFloat(LaneFloat v)
{
	printf("[%f]", v.value);
}

inline void printlnLaneFloat(LaneFloat v)
{
	printf("[%f]", v.value);
}

inline void printLaneVector(LaneVector v)
{
	printf("[%f, %f, %f]", v.x.value, v.y.value, v.z.value);
}

inline void printlnLaneVector(LaneVector v)
{
	printf("[%f, %f, %f]\n", v.x.value, v.y.value, v.z.value);
}*/
