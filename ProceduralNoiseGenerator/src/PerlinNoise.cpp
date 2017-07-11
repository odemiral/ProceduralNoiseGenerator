#include "PerlinNoise.h"




//in 2d there are 8 directions, up, down, left, right, upleft, upright, downleft, downright
const int8_t PerlinNoise::m_grad2DTable[8][2] = {
		{ 0, 1 }, { 1, 1 }, { 1, 0 }, { 1, -1 },
		{ 0, -1 }, { -1, -1 }, { -1, 0 }, { -1, 1 }
};

const int8_t PerlinNoise::m_grad3DTable[16][3] = {
		{ 1, 1, 0 }, { -1, 1, 0 }, { 1, -1, 0 }, { -1, -1, 0 },
		{ 1, 0, 1 }, { -1, 0, 1 }, { 1, 0, -1 }, { -1, 0, -1 },
		{ 0, 1, 1 }, { 0, -1, 1 }, { 0, 1, -1 }, { 0, -1, -1 },
		{ 1, 1, 0 }, { 0, -1, 1 }, { -1, 1, 0 }, { 0, -1, -1 }
};

const int8_t PerlinNoise::m_grad4DTable[32][4] = {
		{ 0, 1, 1, 1 }, { 0, 1, 1, -1 }, { 0, 1, -1, 1 }, { 0, 1, -1, -1 },
		{ 0, -1, 1, 1 }, { 0, -1, 1, -1 }, { 0, -1, -1, 1 }, { 0, -1, -1, -1 },
		{ 1, 0, 1, 1 }, { 1, 0, 1, -1 }, { 1, 0, -1, 1 }, { 1, 0, -1, -1 },
		{ -1, 0, 1, 1 }, { -1, 0, 1, -1 }, { -1, 0, -1, 1 }, { -1, 0, -1, -1 },
		{ 1, 1, 0, 1 }, { 1, 1, 0, -1 }, { 1, -1, 0, 1 }, { 1, -1, 0, -1 },
		{ -1, 1, 0, 1 }, { -1, 1, 0, -1 }, { -1, -1, 0, 1 }, { -1, -1, 0, -1 },
		{ 1, 1, 1, 0 }, { 1, 1, -1, 0 }, { 1, -1, 1, 0 }, { 1, -1, -1, 0 },
		{ -1, 1, 1, 0 }, { -1, 1, -1, 0 }, { -1, -1, 1, 0 }, { -1, -1, -1, 0 }
};

/* Technically it's 256, doubled size to prevent overflowing. if you decide to use anything larger than 512, change table type from u_int8 to int */
//const int PerlinNoise::TABLE_SIZE = 512;
//const int PerlinNoise::TABLE_LIMIT = TABLE_SIZE - 1;



PerlinNoise::PerlinNoise(int octaves, float frequency, float amplitude, int seed, float lacunarity, float gain)
{
	m_seed = seed;
	m_octaves = octaves;
	m_freq = frequency;
	m_amp = amplitude;
	m_lacunarity = lacunarity;
	m_gain = gain;

	srand(m_seed);
	//std::iota(std::begin(m_permTable), std::end(m_permTable), 0); //should init from 0 to 255
	//std::random_shuffle(std::begin(m_permTable), std::end(m_permTable)); //shuffle the permTable

	//for every noise function we'll have a new permutation table.
	std::iota(std::begin(m_permTable), std::end(m_permTable) - TABLE_SIZE / 2, 0); //should init from 0 to 255
	std::random_shuffle(std::begin(m_permTable), std::end(m_permTable) - (TABLE_SIZE / 2)); //shuffle the permTable
	std::copy(std::begin(m_permTable), (std::end(m_permTable) - (TABLE_SIZE / 2)), std::end(m_permTable) - (TABLE_SIZE / 2));
}

/* Fractional Brownian Motion
* given x,y coordinates (usually pixels) return the noise for that pixel.
* In most cases you want to call this function for every pixel in a map.
*  for x to mapx
*        for y to mapy
*           float pixNoise = perlinNoise2D(x,y);
* amplitude and frequency are updated inside this function.
*/
float PerlinNoise::perlinNoise2D(float x, float y) const
{
	float tempAmp = m_amp;
	float tempFreq = m_freq;
	float pixNoise = 0.0f;
	for (int i = 0; i < m_octaves; ++i) {
		pixNoise += noise2D(x * tempFreq, y * tempFreq) * tempAmp;
		tempFreq *= m_lacunarity;
		tempAmp *= m_gain;
	}
	return pixNoise;
}

/* given  x,y,z, generate a perlin noise for that set of coordinates */
float PerlinNoise::perlinNoise3D(float x, float y, float z) const
{
	float tempAmp = m_amp;
	float tempFreq = m_freq;
	float pixNoise = 0.0f;
	for (int i = 0; i < m_octaves; ++i) {
		pixNoise += noise3D(x * tempFreq, y * tempFreq, z * tempFreq) * m_amp;
		tempFreq *= m_lacunarity;
		tempAmp *= m_gain;
	}
	return pixNoise;
}


/* Private Functions */

/*Generates a 2D perlin noise, this function should be called inside perlinNoise2D for each octave. */
float PerlinNoise::noise2D(float x, float y) const
{

	//integer part of x
	int xInt = NoiseHelperFunctions::fastFloor(x);
	int yInt = NoiseHelperFunctions::fastFloor(y);


	float noiseVec[4];
	uint8_t gradVec[4];

	//fraction parts ones that fastFloor and freq 
	float xFraction = x - (float)xInt;
	float yFraction = y - (float)yInt;

	//wrap-around
	xInt &= TABLE_LIMIT;
	yInt &= TABLE_LIMIT;

	/*picking random points from gradient table. */
	//if xInt == 255 
	gradVec[0] = m_permTable[(m_permTable[xInt] + yInt)];
	gradVec[1] = m_permTable[(m_permTable[xInt + 1] + yInt)];
	gradVec[2] = m_permTable[(m_permTable[xInt] + yInt + 1)];
	gradVec[3] = m_permTable[(m_permTable[xInt + 1] + yInt + 1)];

	/* get noises for given points */
	noiseVec[0] = gradDot2D(m_permTable[gradVec[0]], xFraction, yFraction);
	noiseVec[1] = gradDot2D(m_permTable[gradVec[1]], xFraction - 1.0f, yFraction);
	noiseVec[2] = gradDot2D(m_permTable[gradVec[2]], xFraction, yFraction - 1.0f);
	noiseVec[3] = gradDot2D(m_permTable[gradVec[3]], xFraction - 1.0f, yFraction - 1.0f);

	/* fade curve vals */
	float xFade = fade(xFraction);
	float yFade = fade(yFraction);

	//interpolation around x-axis
	const float interx1 = linearInterpolate(noiseVec[0], noiseVec[1], xFade);
	const float interx2 = linearInterpolate(noiseVec[2], noiseVec[3], xFade);

	//interpolation around y-axis
	const float intery = linearInterpolate(interx1, interx2, yFade);
	return intery;
	//return interpolate2d(noiseAA, noiseBA, xFade, noiseAB, noiseBB, yFade);

}

/*Generates a 3D perlin noise, this function should be called inside perlinNoise2D for each octave. */

float PerlinNoise::noise3D(float x, float y, float z) const
{

	/* Gradient table */
	uint8_t gradVec[4];

	/* store noises for points around unit-sphere */
	float noiseVec[8];
	
	//get the integer part of the num, substract 1 if value is negative.
	int xInt = NoiseHelperFunctions::fastFloor(x);
	int yInt = NoiseHelperFunctions::fastFloor(y);
	int zInt = NoiseHelperFunctions::fastFloor(z);

	//fraction parts ones that fastFloor and freq 
	float xFraction = x - (float)xInt;
	float yFraction = y - (float)yInt;
	float zFraction = z - (float)zInt;

	//wrap-around
	xInt &= TABLE_LIMIT;
	yInt &= TABLE_LIMIT;
	zInt &= TABLE_LIMIT;

	/*picking random points from gradient table. Note: 4 corners are enough to get the other 4. */
	gradVec[0] = m_permTable[(m_permTable[xInt] + yInt) & TABLE_LIMIT];
	gradVec[1] = m_permTable[(m_permTable[xInt + 1] + yInt) & TABLE_LIMIT];
	gradVec[2] = m_permTable[(m_permTable[xInt] + yInt + 1) & TABLE_LIMIT];
	gradVec[3] = m_permTable[(m_permTable[xInt + 1] + yInt + 1) & TABLE_LIMIT];


	/* get noises for given points */
	noiseVec[0] = gradDot3D(m_permTable[(gradVec[0] + zInt) & TABLE_LIMIT], xFraction, yFraction, zFraction);
	noiseVec[1] = gradDot3D(m_permTable[(gradVec[1] + zInt) & TABLE_LIMIT], xFraction - 1, yFraction, zFraction);
	noiseVec[2] = gradDot3D(m_permTable[(gradVec[2] + zInt) & TABLE_LIMIT], xFraction, yFraction - 1, zFraction);
	noiseVec[3] = gradDot3D(m_permTable[(gradVec[3] + zInt) & TABLE_LIMIT], xFraction - 1, yFraction - 1, zFraction);

	noiseVec[4] = gradDot3D(m_permTable[(gradVec[0] + zInt + 1) & TABLE_LIMIT], xFraction, yFraction, zFraction - 1);
	noiseVec[5] = gradDot3D(m_permTable[(gradVec[1] + zInt + 1) & TABLE_LIMIT], xFraction - 1, yFraction, zFraction - 1);
	noiseVec[6] = gradDot3D(m_permTable[(gradVec[2] + zInt + 1) & TABLE_LIMIT], xFraction, yFraction - 1, zFraction - 1);
	noiseVec[7] = gradDot3D(m_permTable[(gradVec[3] + zInt + 1) & TABLE_LIMIT], xFraction - 1, yFraction - 1, zFraction - 1);

	/* fade curve vals */
	float xFade = fade(xFraction);
	float yFade = fade(yFraction);
	float zFade = fade(zFraction);


	//interpolate points around unit-circle

	//interpolation around xy-axis
	const float interxy1 = linearInterpolate(noiseVec[0], noiseVec[1], xFade);
	const float interxy2 = linearInterpolate(noiseVec[2], noiseVec[3], xFade);

	//interpolation around xz-axis
	const float interxz1 = linearInterpolate(noiseVec[4], noiseVec[5], xFade);
	const float interxz2 = linearInterpolate(noiseVec[6], noiseVec[7], xFade);

	//interpolation around y-axis
	const float intery1 = linearInterpolate(interxy1, interxy2, yFade);
	const float intery2 = linearInterpolate(interxz1, interxz2, yFade);

	//interpolation around z-axis
	const float interz = linearInterpolate(intery1, intery2, zFade);

	return interz;
}

/* NOT FINISHED YET! */
#ifdef true
float PerlinNoise::noise4D(float x, float y, float z, float t) const
{
	/* Gradient table */
	uint8_t gradVec[4];

	/* store noises for points around unit-3-shpere */
	float noiseVec[16];
}
#endif

/* skew gradient directions away from coordinate axes and long diagonals. */
float PerlinNoise::fade(float t) const
{
	// 6t^5 - 15t^4 + 10t^3  because when t = 0 || t = 1, 1st and 2nd derivatives are 0, reduces the num of artifacts like bump maps.
	//http://mrl.nyu.edu/~perlin/paper445.pdf
	return (t * t * t * (t * (6.0f * t - 15.0f) + 10.0f));
}


/* given hash, x, y values takes the dot product between pseudo-random gradient vector and a point around unit-circle (2d) */
float PerlinNoise::gradDot2D(const uint8_t hash, const float x, const float y) const
{
	const uint8_t newHash = hash & 7;
	return x * m_grad2DTable[newHash][0] + y * m_grad2DTable[newHash][1];
}

/* given hash, x, y, z values takes the dot product between pseudo-random gradient vector and a point around unit-sphere (3d) */
/* TODO: change hash to const uint8_t */
float PerlinNoise::gradDot3D(const uint8_t hash, const float x, const float y, const float z) const
{
	const uint8_t newHash = hash & 15; //use the first 4 bits of the hash value, and use it randomly pick gradient from gradient3d table 
	//(note: random here doesn't mencessarily mean pseudo random, just like with permutation tables it should return the same value for same hash, we want consistency)
	return x * m_grad3DTable[newHash][0] + y * m_grad3DTable[newHash][1] + z * m_grad3DTable[newHash][2];
}

/* given hash, x, y, z, t values takes  dot product between pseudo-random gradient vector and a point around unit-3-sphere (4d)
* Note: 4D usually used to generate animation on procedurally generated 3d content (like ocean, volumetric clouds, smoke simulations)
*/
float PerlinNoise::gradDot4D(const uint8_t hash, const float x, const float y, const float z, const float t) const
{
	const uint8_t newHash = hash & 31;
	return x * m_grad4DTable[newHash][0] + y * m_grad4DTable[newHash][1] + z * m_grad4DTable[newHash][2] + t * m_grad4DTable[newHash][3];
}


/* Interpolate funct takes funct ptr to decide which interpolate funct.
* interpolate2d, 3d, 4d
* interpolate type will be funct ptr later.
*/
float PerlinNoise::interpolate2d(const float leftPoint1, const float rightPoint1, const float valX, const float leftPoint2, const float rightPoint2, const float valY) const
{
	//interpolation around x-axis
	const float x1 = linearInterpolate(leftPoint1, rightPoint1, valX);
	const float x2 = linearInterpolate(leftPoint2, rightPoint2, valX);
	//interpolation around y-axis
	return  linearInterpolate(x1, x2, valY);;
	//return linearInterpolate(linearInterpolate(leftPoint1, rightPoint1, valX), linearInterpolate(leftPoint2, rightPoint2, valX), valY);
}

/* Smoothness of interpolation from smoothest to roughest (also from slowest to fastest)
1) Cubic
2) Cosine
3) Linear
More graphs can be found http://paulbourke.net/miscellaneous/interpolation/
*/
/* interpolation between points @leftPoint and @rightPoint
* @val defines the pos of the estimation on interpolated line. Ranges from 0-1, 0 = val is leftPoint, 1 = val is rightPoint
* Note that these are all 1D interpolate functions,
*/
float PerlinNoise::linearInterpolate(const float leftPoint, const float rightPoint, const float val) const
{
	return leftPoint + val * (rightPoint - leftPoint); //1 less arithmetic op

	//return ((1 - val) * leftPoint + val * rightPoint);
}

float PerlinNoise::cosineInterpolate(const float leftPoint, const float rightPoint, const float val) const
{
	const float ft = val * M_PI;
	const float f = (1 - cos(ft)) * 0.5f; //fastCos produce blocky textures
	return  leftPoint * (1 - f) + rightPoint * f;
}

/* unlike linear and cosine, cubic requires 2 extra points, any point before left point, any point after right point. */
float PerlinNoise::cubicInterpolate(const float beforeLeftPoint, const float leftPoint, const float rightPoint, const float afterRightPoint, const float val) const
{
	const float P = (afterRightPoint - rightPoint) - (beforeLeftPoint - leftPoint);
	const float Q = (beforeLeftPoint - leftPoint) - P;
	const float R = rightPoint - beforeLeftPoint;
	return P * (val * val * val) + Q * (val * val) + R * val + leftPoint; //P*val^3 + Q*val^2 + R*val + leftPoint
}



/* SETTERS */
void PerlinNoise::setOctaves(int octaves)
{
	m_octaves = octaves;
}
void PerlinNoise::setFrequency(float frequency)
{
	m_freq = frequency;
}
void PerlinNoise::setAmplitude(float amplitude)
{
	m_amp = amplitude;
}
void PerlinNoise::setSeed(int seed)
{
	m_seed = seed;
}
void PerlinNoise::setLacunarity(float lacunarity)
{
	m_lacunarity = lacunarity;
}
void PerlinNoise::setGain(float gain)
{
	m_gain = gain;
}


/* GETTERS */
int PerlinNoise::getOctaves() const
{
	return m_octaves;
}
float PerlinNoise::getFrequency() const
{
	return m_freq;
}
float PerlinNoise::getAmplitude() const
{
	return m_amp;
}
int PerlinNoise::getSeed() const
{
	return m_seed;
}
float PerlinNoise::getLacunarity() const
{
	return m_lacunarity;
}
float PerlinNoise::getGain() const
{
	return m_gain;
}



PerlinNoise::~PerlinNoise()
{
}


