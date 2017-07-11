#pragma once
#pragma warning( disable : 4996 ) //for std::copy (visualC++ uses checked iterator for security reasons)
#include <cmath>
#include <iostream> //DELETE
#include <algorithm> //random_shuffle, std::max, and std::sort
#include <stdint.h> //uint8_t
#include <cfloat> // FLT_MAX 


#include "NoiseHelperFunctions.h"
using namespace NoiseHelperFunctions;

/*TODO:
1) Once you're done compare it against std::set to see if it's better to use that (log(N) with overhead ) vs O(Nlog(N))
2) Use variadic functions for dealing with multi dimensional distance functions
*/

/* I believe worley noise is bit more intricate to explain clearly in few sentences, if you like to know how it works, take a look at http://www.rhythmiccanvas.com/research/papers/worley.pdf
 * Worley Noise output will be different based on the distance metric used.
 * If you know any cool metric you know that I haven't implemented, feel free to send me an e-mail.
 * Unlike my other perlin and simplex noise implementations, I decided to use hash function instad of a permutation table.
 */

/*NOTE: Unlike Simplex, or Perlin, I won't be implementing 4D version. It's pretty much the same as 3d but instead of doing 26 checks,
* You need to do 80 checks, since I can't think of any reason for anyone to use 4D on CPU rather GPU, I think it would be a waste of time to implement 4D.
* If you implement one, you're more than welcome to send a merge request.
*/

class WorleyNoise
{
public:
	/* Supported Distance functions, make sure to add any new distance function you implement. */
	enum distanceFunctions{ QUADRATIC, EUCLIDEAN, EUCLIDEAN_SQUARED, MANHATTAN, MINKOWSKI, CITYBLOCK };
	/* Types of noise represented up to F3, these are only meant to help you visualize the effects, you can do whatever you want with feature points. */
	enum noiseTypes{
		FIRST_ORDER, SECOND_ORDER, THIRD_ORDER,
		SUM_F1F2, SUM_F1F2F3, SUM_F2F3, SUM_F1F3,
		DELTA_F2F1, DELTA_F3F2, DELTA_F3F1, DELTA_F3F2F1,
		MULTIPLY_F1F2, MULTIPLY_F2F3, MULTIPLY_F1F3, MULTIPLY_F1F2F3,
		DIVIDE_F2F1, DIVIDE_F3F2, DIVIDE_F3F1, DIVIDE_F3F2F1
	};


	WorleyNoise() = delete;
	WorleyNoise(int octaves, float frequency, float amplitude, int seed, float lacunarity = 2.0f, float gain = 0.65f, distanceFunctions funct = EUCLIDEAN, noiseTypes noiseType = DELTA_F2F1);

	~WorleyNoise();

	/* fractal brownian motion 2D */
	float worleyNoise2D(float x, float y);
	/* fractal brownian motion 3D */
	float worleyNoise3D(float x, float y, float z);


	float turbulanceWorleyNoise2D(float x, float y);
	float turbulanceWorleyNoise3D(float x, float y, float z);

	int getOctaves() const;
	float getFrequency() const;
	float getAmplitude() const;
	int getSeed() const;
	float getLacunarity() const;
	float getGain() const;
	distanceFunctions getDistanceFunction() const;
	noiseTypes getNoiseType() const;


	void setOctaves(int octaves);
	void setFrequency(float frequency);
	void setAmplitude(float amplitude);
	void setSeed(int seed);
	void setLacunarity(float lacunarity);
	void setGain(float gain);
	void setDistanceFunction(distanceFunctions distFunct);
	void setNoiseType(noiseTypes noiseType);
private:

	/* Member Variables */
	float m_amp;
	float m_freq;
	int m_octaves;
	int m_seed;
	float m_lacunarity;
	float m_gain;
	distanceFunctions m_func; /* Distance function. */
	noiseTypes m_noiseType;

	//Constants taken from the sample code from Texturing and Modeling: A Procedural Approach (https://www.amazon.ca/Texturing-Modeling-Procedural-David-Ebert/dp/1558608486)
	
	/* If you decide to have more than 5 Feature points, increase density adjustment.
	* Be careful! increasing density adjustment might cause artifacts.
	*/
	static const int NUM_OF_FEATURES = 2; //MUST BE A POSITIVE INTEGER (DO NECESSARY CHECKS WITHIN CONSTRUCTOR)
	static const float DENSITY_ADJUSTMENT;
	static const float INV_DENSITY_ADJ;
	static const int MINKOWSKI_CONST = 4;
	static const float MINKOWSKI_INV;
	static const float SCALE_CONSTANT;

	float m_featurePts[NUM_OF_FEATURES]; // Value for each feature points
	/* m_delta is the difference between the sample coordinates and  (NUM_OF_FEATURES)th closest feature coordinates.
	* This is faster than dynamically allocate, and resize the array if you ever want to generate multi dimeinsion with same object.
	*/
	float m_delta2D[NUM_OF_FEATURES][2];
	float m_delta3D[NUM_OF_FEATURES][3];

	/* Function ptrs to select the distance funct */
	float (WorleyNoise::*m_distanceFunct2D)(const float(&newPos)[2]);
	float (WorleyNoise::*m_distanceFunct3D)(const float(&newPos)[3]);


	/* Taken from Texture & Modeling: A Procedural Approach 3rd edition.
	* lookup table to determine how many feature points should be in each spatial cube.
	* random indexed num will approximate possion distribution of mean density 2.5
	*/
	static uint8_t m_poissionTable[256];

	/* 2D requires to check 4 vertices, and 4 edges */
	float noise2D(float x, float y);
	/* 3D requires to check 8 vertices, 12 edges, and 6 faces. */
	float noise3D(float x, float y, float z);


	/* Sort square of samples into the current
	* We'll be dealing with the int parts of coordinates.
	*/
	void addSamples2D(int x, int y, const float(&newPos)[2]);
	void addSamples3D(int x, int y, int z, const float(&newPos)[3]);

	//Noise types that are currently supported 
	float noiseType();
	void pickDistanceFunct();
	

	/* 2D  Distance functions */
	inline float quadraticDistance(const float(&deltas)[2])
	{
		return deltas[0] * deltas[0] + deltas[0] * deltas[1] + deltas[1] * deltas[1];
	}
	inline float manhattanDistance(const float(&deltas)[2])
	{
		return fabs(deltas[0]) + fabs(deltas[1]);
	}
	inline float euclideanDistance(const float(&deltas)[2])
	{
		return NoiseHelperFunctions::fastSqrt((deltas[0] * deltas[0]) + (deltas[1] * deltas[1]));
	}
	inline float euclideanDistanceSquared(const float(&deltas)[2])
	{
		return (deltas[0] * deltas[0]) + (deltas[1] * deltas[1]);
	}
	//(abs (deltas[0]) ^ minkowskiConst + abs(deltas[1]) ^ minkowskiConst) ^ (1 / minkowskiConst)
	inline float minkowskiDistance(const float(&deltas)[2])
	{
		return (powf((powf(abs(deltas[0]), MINKOWSKI_CONST) + powf(abs(deltas[1]), MINKOWSKI_CONST)), MINKOWSKI_INV));
	}
	inline float cityBlockDistance(const float(&deltas)[2])
	{
		float dist = std::max(fabs(deltas[0]), fabs(deltas[1]));
		dist *= dist; //may not need this
		return dist;
	}

	/* 3D  Distance functions */
	inline float quadraticDistance(const float(&deltas)[3])
	{
		//x² + y² + z² + xy + xz + yz
		return (deltas[0] * deltas[0]) + (deltas[1] * deltas[1]) + (deltas[2] * deltas[2]) + (deltas[0] * deltas[1]) + (deltas[0] * deltas[2]) + (deltas[1] * deltas[2]);
	}
	inline float manhattanDistance(const float(&deltas)[3])
	{
		return fabs(deltas[0]) + fabs(deltas[1]) + fabs(deltas[2]);
	}
	inline float euclideanDistance(const float(&deltas)[3])
	{
		return NoiseHelperFunctions::fastSqrt((deltas[0] * deltas[0]) + (deltas[1] * deltas[1]) + (deltas[2] * deltas[2]));
	}
	inline float euclideanDistanceSquared(const float(&deltas)[3])
	{
		return (deltas[0] * deltas[0]) + (deltas[1] * deltas[1]) + (deltas[2] * deltas[2]);
	}
	inline float minkowskiDistance(const float(&deltas)[3])
	{
		return (powf((powf(abs(deltas[0]), MINKOWSKI_CONST) + powf(abs(deltas[1]), MINKOWSKI_CONST) + powf(abs(deltas[2]), MINKOWSKI_CONST)), MINKOWSKI_INV));
	}
	inline float cityBlockDistance(const float(&deltas)[3])
	{
		float dist = std::max(std::max(fabs(deltas[0]), fabs(deltas[1])), fabs(deltas[2]));
		dist *= dist;
		return dist;
	}

};

