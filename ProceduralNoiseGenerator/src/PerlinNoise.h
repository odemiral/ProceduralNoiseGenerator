#pragma once

#include <iostream>
#include <algorithm> //random_shuffle
#include <numeric> //iota
#include <stdint.h> //uint8_t

#define _USE_MATH_DEFINES
#include <math.h>

/* TODO: gradientDot functions can be optimized further using SIMD
* TODO: Use XOR operator instead of doubling the size of the table, see how it effects performance. (like http://en.wikipedia.org/wiki/Pearson_hashing)
*  TODO: function ptr for interpolation functions.
* This is a CPU implementation, for any real time rendering, you wanna implement it on GPU.
*
*/
/*
* Octaves: Number of interpolation pass done on the noise, it effectively controls level of detail. Higher the octave, more detailed the image looks.

* Frequency: Determines number of terrain features and size of these features. Decreasing it will decrease number of features and increase size of the features.
There might not be enough pixels to show the details, in fact if frequency is > resolution of the heightmap extra details won't have any affect.

* Amplitude: Amplitude determines the max value for the noise on 3d, amplitude can affect the height of the terrain.
For each octave, amplitude shrinks, when amplitudes are too small, there is no point adding more noise functions since it won't have any significant effect.

* Seed: random seed to pass to srand() which will effect the shuffle pattern (same seed will result in same shuffle pattern), you might want to play around with this value to find interesting patterns.

* Lacunarity: Determines how frequency grows. Usually for each octave, frequency doubles in value.

* Gain: Determines how much the amplitude shrink. Usually for each octave, amplitude halved in value.

* You can see the results of how gain and lacunarity affect signals here: http://freespace.virgin.net/hugo.elias/models/m_perlin.htm
*/

#include "NoiseHelperFunctions.h"

//using namespace NoiseHelperFunctions;
using namespace NoiseHelperFunctions;

class PerlinNoise
{
public:
	PerlinNoise() = delete;
	PerlinNoise(int octaves, float frequency, float amplitude, int seed, float lacunarity = 2.0f, float gain = 0.5f);
	
	~PerlinNoise();

	float perlinNoise2D(float x, float y) const;
	float perlinNoise3D(float x, float y, float z) const;

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


	/*NOTE: All of these private functions were originally inline and it certainly does increase the performance, however for readability reasons I decided to not use inline and thus declare the definition in .cpp file
	 * As I mentioned in my previous comments if you want performance for real time applications, you should opt-in for GPU implementations rathere than this one, but if you still want more performance, then go ahead and use inline
	 */
private:
	float noise2D(float x, float y) const;
	float noise3D(float x, float y, float z) const;
	//const int fastFloor(const float num) const;
	float fade(float t) const; 
	float gradDot2D(const uint8_t hash, const float x, const float y) const;
	float gradDot3D(const uint8_t hash, const float x, const float y, const float z) const;
	float gradDot4D(const uint8_t hash, const float x, const float y, const float z, const float t) const;
	float interpolate2d(const float leftPoint1, const float rightPoint1, const float valX, const float leftPoint2, const float rightPoint2, const float valY) const;
	float linearInterpolate(const float leftPoint, const float rightPoint, const float val) const;
	float cosineInterpolate(const float leftPoint, const float rightPoint, const float val) const;
	float cubicInterpolate(const float beforeLeftPoint, const float leftPoint, const float rightPoint, const float afterRightPoint, const float val) const;


	/* Technically it's 256, doubled size to prevent overflowing. if you decide to use anything larger than 512, change table type from u_int8 to int */
	const static int TABLE_SIZE = 512;
	static const int TABLE_LIMIT = TABLE_SIZE - 1;

	/* Member variables */
	float m_amp;
	float m_freq;
	int m_octaves;
	int m_seed;
	uint8_t m_permTable[TABLE_SIZE];
	float m_lacunarity;
	float m_gain;

	/* Gradient tables */
	static const int8_t m_grad2DTable[8][2];
	static const int8_t m_grad3DTable[16][3];
	static const int8_t m_grad4DTable[32][4];

};
