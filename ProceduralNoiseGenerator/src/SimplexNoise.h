#pragma once
#include <iostream>
#include <algorithm> //random_shuffle
#include <numeric> //iota
#include <stdint.h> //uint8_t

#define _USE_MATH_DEFINES
#include <math.h>

#include "NoiseHelperFunctions.h"

/* Optimizations are partially based on http://webstaff.itn.liu.se/~stegu/simplexnoise/SimplexNoise.java
* Further optimization can be done by removing 2d table and using x,y of 3d for 2d instead.
* You can also try to store temp variables in noise functions to reduce allocation overhead.
* TODO: I think using & operation every time on array is less efficient than doubling the size of the array at construction time.
*/

class SimplexNoise
{
public:
	SimplexNoise() = delete;
	SimplexNoise(int octaves, float frequency, float amplitude, int seed, float lacunarity = 2.0f, float gain = 0.65f);

	~SimplexNoise();

	float simplexNoise2D(float x, float y) const;
	float simplexNoise3D(float x, float y, float z) const;
	float simplexNoise4D(float x, float y, float z, float t) const;

	int getOctaves() const;
	float getFrequency() const;
	float getAmplitude() const;
	int getSeed() const;
	float getLacunarity() const;
	float getGain() const;

	void setOctaves(int octaves);
	void setFrequency(float frequency);
	void setAmplitude(float amplitude);
	void setSeed(int seed);
	void setLacunarity(float lacunarity);
	void setGain(float gain);

private:

	/* Member variables */
	float m_amp;
	float m_freq;
	int m_octaves;
	int m_seed;
	float m_lacunarity;
	float m_gain;

	/* Technically it's 256, doubled size to prevent overflowing. if you decide to use anything larger than 512, change table type from u_int8 to int */
	static const int TABLE_SIZE = 512;
	static const int TABLE_LIMIT = TABLE_SIZE - 1;


	/* Skewing factors for 2,3,4 dimensions (skewing factor for Nth dimension = (sqrt(N+1) - 1) / N
	* Needed for squishing n-dimensional hypercubic lattices.
	*/
	static const float skewF2D;
	static const float skewF3D;
	static const float skewF4D;

	/* Unskewing factors for 2,3,4 dimensions (unskewing factor for Nth dimension = N + 1 - sqrt(N+1) /(N*(N+1)) */
	static const float unskewF2D;
	static const float unskewF3D;
	static const float unskewF4D;


	uint8_t m_permTable[TABLE_SIZE];

	//static lookup tables
	static const float m_grad2DTable[8][2];
	static const uint8_t m_grad3DTable[16][3];
	static const float m_grad4DTable[32][4];
	static const float m_simplexTable[64][4];

	//single simplex noise used by 
	float noise2D(float x, float y) const;
	float noise3D(float x, float y, float z) const;
	float noise4D(const float x, const float y, const float z, const float t) const;

	//helper functions that finds the dot product between 2 gradient vectors.
	float gradDot2D(const uint8_t hash, const float x, const float y) const;
	float gradDot3D(const uint8_t hash, const float x, const float y, const float z) const;
	float gradDot4D(const uint8_t hash, const float x, const float y, const float z, const float t) const;

};

