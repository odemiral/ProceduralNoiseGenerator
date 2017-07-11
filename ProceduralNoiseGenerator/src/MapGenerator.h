/* Simple 2D Map Generator using various procedural noises
 * I mainly used this to debug my code, Procedural noise bugs can be quite hard to track without a proper visualization
*/

#pragma warning( disable : 4996 )
#pragma once

#include <iostream>
#include <string>
#include <fstream>

#define _USE_MATH_DEFINES
#include <math.h>


#include "PerlinNoise.h"
#include "SimplexNoise.h"
#include "WorleyNoise.h"


using namespace std;

/* RGB color structure */
struct color{
	//v[0]=red, v[1]=green, v[2]=blue
	unsigned char v[3];

	color(unsigned char r, unsigned char g, unsigned char b){
		v[0] = r;
		v[1] = g;
		v[2] = b;
	}
	friend color operator- (color const& lhs, color const& rhs) {
		unsigned char r = lhs.v[0] - rhs.v[0];
		unsigned char g = lhs.v[1] - rhs.v[1];
		unsigned char b = lhs.v[2] - rhs.v[2];
		return color(r, g, b);
	}
	friend color operator* (color const& lhs, color const& rhs) {
		unsigned char r = lhs.v[0] * rhs.v[0];
		unsigned char g = lhs.v[1] * rhs.v[1];
		unsigned char b = lhs.v[2] * rhs.v[2];
		return color(r, g, b);
	}
	friend color operator* (float lhs, color const& rhs) {
		unsigned char r = lhs * rhs.v[0];
		unsigned char g = lhs * rhs.v[1];
		unsigned char b = lhs * rhs.v[2];
		return color(r, g, b);
	}

	friend color operator+ (color const& lhs, color const& rhs) {
		unsigned char r = lhs.v[0] + rhs.v[0];
		unsigned char g = lhs.v[1] + rhs.v[1];
		unsigned char b = lhs.v[2] + rhs.v[2];
		return color(r, g, b);
	}
};

/* MAP GENERATOR Functions Overview
* Filler Functions: They fill the map with specified procedural noises
*/


/* TODO: Convert m_map to vector and make it 3d, for efficiency you can make your own MapContainer Class to deal with the data.
*/
class MapGenerator
{
public:
	enum Noise{ PERLIN, SIMPLEX, WORLEY };
	MapGenerator() = delete;
	MapGenerator(float hgrid, float vgrid, string oname);
	~MapGenerator();
	//interpolation functions

	void createPerlinNoise(int octaves, float frequency, float amplitude, int seed, float lacunarity = 2.0f, float gain = 0.5f);
	void createSimplexNoise(int octaves, float frequency, float amplitude, int seed, float lacunarity = 2.0f, float gain = 0.5f);
	void createWorleyNoise(int octaves, float frequency, float amplitude, int seed, float lacunarity = 2.0f, float gain = 0.5f,
		WorleyNoise::distanceFunctions func = WorleyNoise::distanceFunctions::EUCLIDEAN, WorleyNoise::noiseTypes noiseType = WorleyNoise::noiseTypes::DELTA_F2F1);

	//given noise type, fills the map with 2-D of that noise.
	void fillMapWithNoise2D(Noise noise);


	//clears data inside the map and fills it with 0s
	void flushMap();

	void printNoiseInfo() const;

	void printColoredImage(string imageName) const;
	void printGrayscaleImage(string imageName) const;
private:
	/* Properties used by the Map */
	int m_hgrid, m_vgrid;
	float m_min, m_max;  //min and max noise values in the map, depending on the application it could be used as an elevation data or to determine the coloring differences of the map.
	string m_mapName;
	float **m_map;

	WorleyNoise *m_wnoise;
	//SimplexNoise *m_snoise;
	PerlinNoise *m_pnoise;
	SimplexNoise *m_snoise;


	inline color lerp(color leftPoint, color rightPoint, float val) const
	{
		const float ft = val * M_PI;
		const float f2 = (1 - cos(ft)) * 0.5f;
		color c2 = (1 - f2) * leftPoint + f2 * rightPoint;
		return c2;
	}
	inline color lerp2(color c1, color c2, float value) const
	{
		return ((1 - value) * c1 + value * c2);
	}

	//Given value x in [A,B] range, scales it to [C,D] range. Useful when you want to convert min,max value of a noise to rgb [0,255] range.
	/* (C/(B-A))(B - x) + (D/(B-A)) (x - A) == C + (D?C) * ((x?A) / (B?A))
	*/
	inline int scale(float x, float A, float B, float C, float D) const
	{
	//return C * (1 - ((x - A) / (B - A))) + D * ((x - A) / (B - A));
	return (int)(C + ((D - C) * ((x - A) / (B - A)))); // simplified
	}
};

