
#pragma once
namespace NoiseHelperFunctions
{
	/* Floors the value to get the surrounding points along unit-circle/sphere/3-sphere */
	inline int fastFloor(const float num)
	{
		return (num >= 0.0f ? (int)num : (int)num - 1);
	}

	/* Donald Knuth's LCG */
	inline void LCG(unsigned long &seed) { seed = 1402024253 * seed + 586950981; }

	/* fast but not-very accurate way of calculating sqrt */
	inline float fastSqrt(float x)
	{
		unsigned int i = *(unsigned int*)&x;
		// adjust bias
		i += 127 << 23;
		// approximation of square root
		i >>= 1;
		return *(float*)&i;
	}
}