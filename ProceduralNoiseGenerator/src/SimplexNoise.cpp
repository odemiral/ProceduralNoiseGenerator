#include "SimplexNoise.h"


/* CONSTANT VARIABLE INITIALIZATIONS */
//in 2d there are 8 directions, up, down, left, right, upleft, upright, downleft, downright
const float SimplexNoise::m_grad2DTable[8][2] = {
		{ 0, 1 }, { 1, 1 }, { 1, 0 }, { 1, -1 },
		{ 0, -1 }, { -1, -1 }, { -1, 0 }, { -1, 1 }
};

const uint8_t SimplexNoise::m_grad3DTable[16][3] = {
		{ 1, 1, 0 }, { -1, 1, 0 }, { 1, -1, 0 }, { -1, -1, 0 },
		{ 1, 0, 1 }, { -1, 0, 1 }, { 1, 0, -1 }, { -1, 0, -1 },
		{ 0, 1, 1 }, { 0, -1, 1 }, { 0, 1, -1 }, { 0, -1, -1 },
		{ 1, 1, 0 }, { 0, -1, 1 }, { -1, 1, 0 }, { 0, -1, -1 }
};

const float SimplexNoise::m_grad4DTable[32][4] = {
		{ 0, 1, 1, 1 }, { 0, 1, 1, -1 }, { 0, 1, -1, 1 }, { 0, 1, -1, -1 },
		{ 0, -1, 1, 1 }, { 0, -1, 1, -1 }, { 0, -1, -1, 1 }, { 0, -1, -1, -1 },
		{ 1, 0, 1, 1 }, { 1, 0, 1, -1 }, { 1, 0, -1, 1 }, { 1, 0, -1, -1 },
		{ -1, 0, 1, 1 }, { -1, 0, 1, -1 }, { -1, 0, -1, 1 }, { -1, 0, -1, -1 },
		{ 1, 1, 0, 1 }, { 1, 1, 0, -1 }, { 1, -1, 0, 1 }, { 1, -1, 0, -1 },
		{ -1, 1, 0, 1 }, { -1, 1, 0, -1 }, { -1, -1, 0, 1 }, { -1, -1, 0, -1 },
		{ 1, 1, 1, 0 }, { 1, 1, -1, 0 }, { 1, -1, 1, 0 }, { 1, -1, -1, 0 },
		{ -1, 1, 1, 0 }, { -1, 1, -1, 0 }, { -1, -1, 1, 0 }, { -1, -1, -1, 0 }
};

/* Lookup tabke for traversing the simplex around any given point in 4D */
const float SimplexNoise::m_simplexTable[64][4] = {
		{ 0, 1, 2, 3 }, { 0, 1, 3, 2 }, { 0, 0, 0, 0 }, { 0, 2, 3, 1 }, { 0, 0, 0, 0 }, { 0, 0, 0, 0 }, { 0, 0, 0, 0 }, { 1, 2, 3, 0 },
		{ 0, 2, 1, 3 }, { 0, 0, 0, 0 }, { 0, 3, 1, 2 }, { 0, 3, 2, 1 }, { 0, 0, 0, 0 }, { 0, 0, 0, 0 }, { 0, 0, 0, 0 }, { 1, 3, 2, 0 },
		{ 0, 0, 0, 0 }, { 0, 0, 0, 0 }, { 0, 0, 0, 0 }, { 0, 0, 0, 0 }, { 0, 0, 0, 0 }, { 0, 0, 0, 0 }, { 0, 0, 0, 0 }, { 0, 0, 0, 0 },
		{ 1, 2, 0, 3 }, { 0, 0, 0, 0 }, { 1, 3, 0, 2 }, { 0, 0, 0, 0 }, { 0, 0, 0, 0 }, { 0, 0, 0, 0 }, { 2, 3, 0, 1 }, { 2, 3, 1, 0 },
		{ 1, 0, 2, 3 }, { 1, 0, 3, 2 }, { 0, 0, 0, 0 }, { 0, 0, 0, 0 }, { 0, 0, 0, 0 }, { 2, 0, 3, 1 }, { 0, 0, 0, 0 }, { 2, 1, 3, 0 },
		{ 0, 0, 0, 0 }, { 0, 0, 0, 0 }, { 0, 0, 0, 0 }, { 0, 0, 0, 0 }, { 0, 0, 0, 0 }, { 0, 0, 0, 0 }, { 0, 0, 0, 0 }, { 0, 0, 0, 0 },
		{ 2, 0, 1, 3 }, { 0, 0, 0, 0 }, { 0, 0, 0, 0 }, { 0, 0, 0, 0 }, { 3, 0, 1, 2 }, { 3, 0, 2, 1 }, { 0, 0, 0, 0 }, { 3, 1, 2, 0 },
		{ 2, 1, 0, 3 }, { 0, 0, 0, 0 }, { 0, 0, 0, 0 }, { 0, 0, 0, 0 }, { 3, 1, 0, 2 }, { 0, 0, 0, 0 }, { 3, 2, 0, 1 }, { 3, 2, 1, 0 }
};


const float SimplexNoise::skewF2D = 0.36602540378444f;
const float SimplexNoise::skewF3D = 0.33333333333334f;
const float SimplexNoise::skewF4D = 0.30901699437495f;

const float SimplexNoise::unskewF2D = 0.211324865405f;
const float SimplexNoise::unskewF3D = 0.166666666667f;
const float SimplexNoise::unskewF4D = 0.138196601125f;


SimplexNoise::SimplexNoise(int octaves, float frequency, float amplitude, int seed, float lacunarity, float gain)
{
	m_seed = seed;
	m_octaves = octaves;
	m_freq = frequency;
	m_amp = amplitude;
	m_lacunarity = lacunarity;
	m_gain = gain;

	srand(m_seed);

	std::iota(std::begin(m_permTable), std::end(m_permTable) - TABLE_SIZE / 2, 0); //should init from 0 to TABLE_SIZE / 2
	std::random_shuffle(std::begin(m_permTable), std::end(m_permTable) - (TABLE_SIZE / 2)); //shuffle the first half
	/* Copy the first half of the array to second half. */
	std::copy(std::begin(m_permTable), (std::end(m_permTable) - (TABLE_SIZE / 2)), std::end(m_permTable) - (TABLE_SIZE / 2));

	//for (int i = 0; i < 256; ++i) { std::cout << int(m_permTable[i]) << ", "; } std::cout << std::endl; system("pause"); std::cout << "__________\n"; for (int i = 256; i < 512; ++i) { std::cout << int(m_permTable[i]) << ", "; } system("pause");
}

float SimplexNoise::simplexNoise2D(float x, float y) const
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

float SimplexNoise::simplexNoise3D(float x, float y, float z) const
{
	float tempAmp = m_amp;
	float tempFreq = m_freq;
	float pixNoise = 0.0f;
	for (int i = 0; i < m_octaves; ++i) {
		pixNoise += noise3D(x * tempFreq, y * tempFreq, z * tempFreq) * tempAmp;
		tempFreq *= m_lacunarity;
		tempAmp *= m_gain;
	}

	return pixNoise;
}

float SimplexNoise::simplexNoise4D(float x, float y, float z, float t) const
{
	float tempAmp = m_amp;
	float tempFreq = m_freq;
	float pixNoise = 0.0f;
	for (int i = 0; i < m_octaves; ++i) {
		pixNoise += noise4D(x * tempFreq, y * tempFreq, z * tempFreq, t * tempFreq) * tempAmp;
		tempFreq *= m_lacunarity;
		tempAmp *= m_gain;
	}
	return pixNoise;
}

float SimplexNoise::noise2D(float x, float y) const
{
	/* store contributions around 3 corners using attenuation function */
	float bContribution, mContribution, tContribution;

	/* Simplex corners for x,y */
	int bCornerX, bCornerY,
		mCornerX, mCornerY,
		tCornerX, tCornerY;

	/* simplex noises for each corner */
	float noiseVec[3] = { 0.0f, 0.0f, 0.0f };

	/* gradiens for each corner */
	uint8_t gradVec[3];

	/* distance values for x,y from each corner */
	float distanceFromBottom[2];
	float distanceFromMid[2];
	float distanceFromTop[2];

	float skew = (x + y) * skewF2D;
	//get bottom and top corners
	bCornerX = NoiseHelperFunctions::fastFloor(x + skew); //bottom corner
	bCornerY = NoiseHelperFunctions::fastFloor(y + skew); //bottom corner
	tCornerX = 1 + bCornerX;
	tCornerY = 1 + bCornerY;

	//compute skew for 2d
	float unskew = (bCornerX + bCornerY) * unskewF2D;
	distanceFromBottom[0] = x - (bCornerX - unskew);
	distanceFromBottom[1] = y - (bCornerY - unskew);

	/*	if distbot x,y is (0,0),(1,1) then midCorner is (0,1)
	else midCorner is (1,0)
	This way there is no need for branches.
	*/
	mCornerX = distanceFromBottom[0] > distanceFromBottom[1];
	mCornerY = distanceFromBottom[0] <= distanceFromBottom[1];

	distanceFromMid[0] = distanceFromBottom[0] - mCornerX + unskewF2D; //float(mCornerX - bCornerX) + unskewF2D;
	distanceFromMid[1] = distanceFromBottom[1] - mCornerY + unskewF2D; // float(mCornerY - bCornerY) + unskewF2D;

	distanceFromTop[0] = distanceFromBottom[0] - 1.0f + unskewF2D + unskewF2D;
	distanceFromTop[1] = distanceFromBottom[1] - 1.0f + unskewF2D + unskewF2D;

	//get the gradients indices
	bCornerX &= TABLE_LIMIT;
	bCornerY &= TABLE_LIMIT;

	gradVec[0] = m_permTable[m_permTable[bCornerY] + bCornerX];
	gradVec[1] = m_permTable[m_permTable[bCornerY + mCornerY] + mCornerX + bCornerX];
	gradVec[2] = m_permTable[m_permTable[tCornerY] + tCornerX];


	/* calculate contributions  (atten coef is 0.5) */
	bContribution = 0.5f - distanceFromBottom[0] * distanceFromBottom[0] - distanceFromBottom[1] * distanceFromBottom[1];
	mContribution = 0.5f - distanceFromMid[0] * distanceFromMid[0] - distanceFromMid[1] * distanceFromMid[1];
	tContribution = 0.5f - distanceFromTop[0] * distanceFromTop[0] - distanceFromTop[1] * distanceFromTop[1];


	//for a given corner if the noise contribution is > 0, apply the contribution to the noise vec. 
	if (bContribution > 0) {
		noiseVec[0] = powf(bContribution, 4.0f) * gradDot2D(gradVec[0], distanceFromBottom[0], distanceFromBottom[1]); //bContribution * bContribution * bContribution * bContribution * 
	}
	if (mContribution > 0) {
		noiseVec[1] = powf(mContribution, 4.0f) * gradDot2D(gradVec[1], distanceFromMid[0], distanceFromMid[1]); //bContribution * bContribution * bContribution * bContribution * 
	}
	if (tContribution > 0) {
		noiseVec[2] = powf(tContribution, 4.0f) * gradDot2D(gradVec[2], distanceFromTop[0], distanceFromTop[1]); //bContribution * bContribution * bContribution * bContribution * 
	}

	return 70.0f * (noiseVec[0] + noiseVec[1] + noiseVec[2]);
}

float SimplexNoise::noise3D(float x, float y, float z) const
{
	/* Simplex noise for each corner */
	float noiseVec[4] = { 0.0f, 0.0f, 0.0f, 0.0f };

	/* Noise contributions around 4 corners of simplex */
	float noiseContrVec[4];

	/* gradients of 4 simplex corners */
	uint8_t gradVec[4];

	/* Vectors storing x,y,z vals for each corner (it should look like a tetrahedron
	No need to store forth corner (last) since it's going to be first corner indices + 1
	Corner values will only be 0 or 1
	*/
	uint8_t firstCornerVec[3], secondCornerVec[3], thirdCornerVec[3]; //TEST IF OVERFLOWS
	/* Each corner vector stores x,y,z values for a given simplex corners */
	float distFromFirstCorner[3], distFromSecondCorner[3], distFromThirdCorner[3], distFromForthCorner[3];

	float skew = (x + y + z) * skewF3D;
	firstCornerVec[0] = NoiseHelperFunctions::fastFloor(x + skew); //skewed x
	firstCornerVec[1] = NoiseHelperFunctions::fastFloor(y + skew); //skewed y
	firstCornerVec[2] = NoiseHelperFunctions::fastFloor(z + skew); //skewed z


	float unskew = (firstCornerVec[0] + firstCornerVec[1] + firstCornerVec[2]) * unskewF3D;
	distFromFirstCorner[0] = x - (firstCornerVec[0] - unskew);
	distFromFirstCorner[1] = y - (firstCornerVec[1] - unskew);
	distFromFirstCorner[2] = z - (firstCornerVec[2] - unskew);



	// THIS PART SEEMS LIKE IT CAN BE EXPLOITED BY USING SOME CLEVER BIT SHIFTING 
	if (distFromFirstCorner[0] < distFromFirstCorner[1]) {
		if (distFromFirstCorner[1] < distFromFirstCorner[2]) {
			secondCornerVec[0] = 0; secondCornerVec[1] = 0; secondCornerVec[2] = 1;
			thirdCornerVec[0] = 0;	thirdCornerVec[1] = 1;	thirdCornerVec[2] = 1;
		}
		else if (distFromFirstCorner[0] < distFromFirstCorner[2]) {
			secondCornerVec[0] = 0; secondCornerVec[1] = 1; secondCornerVec[2] = 0;
			thirdCornerVec[0] = 0;	thirdCornerVec[1] = 1;	thirdCornerVec[2] = 1;
		}
		else {
			secondCornerVec[0] = 0; secondCornerVec[1] = 1; secondCornerVec[2] = 0;
			thirdCornerVec[0] = 1;	thirdCornerVec[1] = 1;	thirdCornerVec[2] = 0;
		}
	}
	else {
		if (distFromFirstCorner[1] >= distFromFirstCorner[2]) {
			secondCornerVec[0] = 1; secondCornerVec[1] = 0; secondCornerVec[2] = 0;
			thirdCornerVec[0] = 1;	thirdCornerVec[1] = 1;	thirdCornerVec[2] = 0;
		}
		else if (distFromFirstCorner[0] >= distFromFirstCorner[2]) {
			secondCornerVec[0] = 1; secondCornerVec[1] = 0; secondCornerVec[2] = 0;
			thirdCornerVec[0] = 1;	thirdCornerVec[1] = 0;	thirdCornerVec[2] = 1;
		}
		else {
			secondCornerVec[0] = 0; secondCornerVec[1] = 0; secondCornerVec[2] = 1;
			thirdCornerVec[0] = 1;	thirdCornerVec[1] = 0;	thirdCornerVec[2] = 1;
		}
	}

	/* Calculate the distance value of x,y,z from 2nd, 3rd and 4th corners */
	distFromSecondCorner[0] = distFromFirstCorner[0] - (secondCornerVec[0] - unskewF3D);
	distFromSecondCorner[1] = distFromFirstCorner[1] - (secondCornerVec[1] - unskewF3D);
	distFromSecondCorner[2] = distFromFirstCorner[2] - (secondCornerVec[2] - unskewF3D);

	distFromThirdCorner[0] = distFromFirstCorner[0] - (thirdCornerVec[0] - 2.0f * unskewF3D);
	distFromThirdCorner[1] = distFromFirstCorner[1] - (thirdCornerVec[1] - 2.0f * unskewF3D);
	distFromThirdCorner[2] = distFromFirstCorner[2] - (thirdCornerVec[2] - 2.0f * unskewF3D);

	distFromForthCorner[0] = distFromFirstCorner[0] - (1.0f - 3.0f * unskewF3D);
	distFromForthCorner[1] = distFromFirstCorner[1] - (1.0f - 3.0f * unskewF3D);
	distFromForthCorner[2] = distFromFirstCorner[2] - (1.0f - 3.0f * unskewF3D);

	/* get the gradients indices */
	firstCornerVec[0] &= TABLE_LIMIT;
	firstCornerVec[1] &= TABLE_LIMIT;
	firstCornerVec[2] &= TABLE_LIMIT;

	gradVec[0] = m_permTable[m_permTable[firstCornerVec[1] + m_permTable[firstCornerVec[2]]] + firstCornerVec[0]];
	gradVec[1] = m_permTable[m_permTable[firstCornerVec[1] + secondCornerVec[1] + m_permTable[firstCornerVec[2] + secondCornerVec[2]]] + firstCornerVec[0] + secondCornerVec[0]];
	gradVec[2] = m_permTable[m_permTable[firstCornerVec[1] + thirdCornerVec[1] + m_permTable[firstCornerVec[2] + thirdCornerVec[2]]] + firstCornerVec[0] + thirdCornerVec[0]];
	gradVec[3] = m_permTable[m_permTable[firstCornerVec[1] + 1 + m_permTable[firstCornerVec[2] + 1]] + firstCornerVec[0] + 1];

	noiseContrVec[0] = 0.5f - (distFromFirstCorner[0] * distFromFirstCorner[0]) - (distFromFirstCorner[1] * distFromFirstCorner[1]) - (distFromFirstCorner[2] * distFromFirstCorner[2]);
	noiseContrVec[1] = 0.5f - (distFromSecondCorner[0] * distFromSecondCorner[0]) - (distFromSecondCorner[1] * distFromSecondCorner[1]) - (distFromSecondCorner[2] * distFromSecondCorner[2]);
	noiseContrVec[2] = 0.5f - (distFromThirdCorner[0] * distFromThirdCorner[0]) - (distFromThirdCorner[1] * distFromThirdCorner[1]) - (distFromThirdCorner[2] * distFromThirdCorner[2]);
	noiseContrVec[3] = 0.5f - (distFromForthCorner[0] * distFromForthCorner[0]) - (distFromForthCorner[1] * distFromForthCorner[1]) - (distFromForthCorner[2] * distFromForthCorner[2]);

	//for a given corner if the noise contribution is > 0, apply the contribution to the noise vec. 
	if (noiseContrVec[0] > 0) {
		noiseVec[0] = pow(noiseContrVec[0], 4.0f) * gradDot3D(gradVec[0], distFromFirstCorner[0], distFromFirstCorner[1], distFromFirstCorner[2]); //bContribution * bContribution * bContribution * bContribution * 
	}
	if (noiseContrVec[1] > 0) {
		noiseVec[1] = pow(noiseContrVec[1], 4.0f) * gradDot3D(gradVec[1], distFromSecondCorner[0], distFromSecondCorner[1], distFromSecondCorner[2]); //bContribution * bContribution * bContribution * bContribution * 
	}
	if (noiseContrVec[2] > 0) {
		noiseVec[2] = pow(noiseContrVec[2], 4.0f) * gradDot3D(gradVec[2], distFromThirdCorner[0], distFromThirdCorner[1], distFromThirdCorner[2]); //bContribution * bContribution * bContribution * bContribution * 
	}
	if (noiseContrVec[3] > 0) {
		noiseVec[3] = pow(noiseContrVec[3], 4.0f) * gradDot3D(gradVec[3], distFromForthCorner[0], distFromForthCorner[1], distFromForthCorner[2]); //bContribution * bContribution * bContribution * bContribution * 
	}

	return 32.0f * (noiseVec[0] + noiseVec[1] + noiseVec[2] + noiseVec[3]);
}


/*This one uses simplex look-up tables to increase performance.  */
float SimplexNoise::noise4D(const float x, const float y, const float z, const float t) const
{
	float noiseVec[5] = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };

	/* Noise contributions around 5 corners of simplex */
	float noiseContrVec[5];

	/* gradients of 5 simplex corners */
	uint8_t gradVec[5];

	/* Vectors storing x,y,z,t vals for each corner (I have no idea what it would look like...)
	No need to store fifth corner (last) since it's going to be all 1s (1,1,1,1)
	*/
	uint8_t firstCornerVec[4], secondCornerVec[4], thirdCornerVec[4], forthCornerVec[4]; //TEST IF OVERFLOWS
	/* Each corner vector stores x,y,z,t values for a given simplex corners */
	float distFromFirstCorner[4], distFromSecondCorner[4], distFromThirdCorner[4], distFromForthCorner[4], distFromFifthCorner[4];

	float skew = (x + y + z + t) * skewF4D;
	firstCornerVec[0] = NoiseHelperFunctions::fastFloor(x + skew); //skewed x
	firstCornerVec[1] = NoiseHelperFunctions::fastFloor(y + skew); //skewed y
	firstCornerVec[2] = NoiseHelperFunctions::fastFloor(z + skew); //skewed z
	firstCornerVec[3] = NoiseHelperFunctions::fastFloor(t + skew); //skewed t

	float unskew = (firstCornerVec[0] + firstCornerVec[1] + firstCornerVec[2] + firstCornerVec[3]) * unskewF4D;
	distFromFirstCorner[0] = x - (firstCornerVec[0] - unskew);
	distFromFirstCorner[1] = y - (firstCornerVec[1] - unskew);
	distFromFirstCorner[2] = z - (firstCornerVec[2] - unskew);
	distFromFirstCorner[3] = t - (firstCornerVec[3] - unskew);

	/* Wraparound */
	firstCornerVec[0] &= TABLE_LIMIT;
	firstCornerVec[1] &= TABLE_LIMIT;
	firstCornerVec[2] &= TABLE_LIMIT;
	firstCornerVec[3] &= TABLE_LIMIT;

	/* given coordinates there are 24 possible simplices, 6 pairs for each 4 coordinates, usually you need to do 6 comparison to find
	the traversal order between simplices, compare the distances between simplices to find the correct order (x,y,z,t) and assign values depending
	value comparison between coordinates.
	if distance from x > dist from y:
	the pair's coordinates has the largest magnitude
	else pair's coordinates = 0
	if dist from x is less than dist from z:
	the pair's coordinates has the 2nd largest mag
	else pair's coordinates = 0
	do that for all pairs and sum them to get the index for simplex lookup table.
	Finally sum them all which will be used as an index for simplex lookups
	Magnitudes are pow of 2 since the size of simplex lookup table is 64 (32 + 16 + 8 + 4 + 2 + 1)

	Update: There is no need to use branches, since values will be between pow2 or 0 we can just sum bool conditions.
	*/
	uint8_t simplexIndex =
		(distFromFirstCorner[0] > distFromFirstCorner[1]) * 32 +
		(distFromFirstCorner[0] > distFromFirstCorner[2]) * 16 +
		(distFromFirstCorner[1] > distFromFirstCorner[2]) * 8 +
		(distFromFirstCorner[0] > distFromFirstCorner[3]) * 4 +
		(distFromFirstCorner[1] > distFromFirstCorner[4]) * 2 +
		(distFromFirstCorner[2] > distFromFirstCorner[4]);

	/* TRY
	uint8_t simplexIndex =
	(distFromFirstCorner[0] > distFromFirstCorner[1]) * 32 |
	(distFromFirstCorner[0] > distFromFirstCorner[2]) * 16 |
	(distFromFirstCorner[1] > distFromFirstCorner[2]) * 8 |
	(distFromFirstCorner[0] > distFromFirstCorner[3]) * 4 |
	(distFromFirstCorner[1] > distFromFirstCorner[4]) * 2 |
	(distFromFirstCorner[2] > distFromFirstCorner[4]);
	*/

	//get 2nd, 3rd, 4th corner coordinates
	/* among 64 entries, only 24 will be possible to reach,
	so we check the simplex table to find possible coordinates (one that's within the reachable 24) for each corner.
	*/
	//for the 2nd corner we check the largest coordinates in simplex table (table vals: 0,1,2,3)
	secondCornerVec[0] = m_simplexTable[simplexIndex][0] >= 3;
	secondCornerVec[1] = m_simplexTable[simplexIndex][1] >= 3;
	secondCornerVec[2] = m_simplexTable[simplexIndex][2] >= 3;
	secondCornerVec[3] = m_simplexTable[simplexIndex][3] >= 3;

	//for the 3rd corner we check the 2nd largest coordinates in the table
	thirdCornerVec[0] = m_simplexTable[simplexIndex][0] >= 2;
	thirdCornerVec[1] = m_simplexTable[simplexIndex][1] >= 2;
	thirdCornerVec[2] = m_simplexTable[simplexIndex][2] >= 2;
	thirdCornerVec[3] = m_simplexTable[simplexIndex][3] >= 2;

	//for the 4th corner we check the 3rd largest coordinates in the table
	forthCornerVec[0] = m_simplexTable[simplexIndex][0] >= 1;
	forthCornerVec[1] = m_simplexTable[simplexIndex][1] >= 1;
	forthCornerVec[2] = m_simplexTable[simplexIndex][2] >= 1;
	forthCornerVec[3] = m_simplexTable[simplexIndex][3] >= 1;
	//since fifth one coordaintes will be 1, simplex should be 0, no need to do any checking.

	/* Calculate distances for 2nd, 3rd, 4th and 5th corners */
	distFromSecondCorner[0] = distFromFirstCorner[0] - (secondCornerVec[0] - unskewF3D);
	distFromSecondCorner[1] = distFromFirstCorner[1] - (secondCornerVec[1] - unskewF3D);
	distFromSecondCorner[2] = distFromFirstCorner[2] - (secondCornerVec[2] - unskewF3D);
	distFromSecondCorner[3] = distFromFirstCorner[3] - (secondCornerVec[3] - unskewF3D);

	distFromThirdCorner[0] = distFromFirstCorner[0] - (thirdCornerVec[0] - 2.0f * unskewF3D);
	distFromThirdCorner[1] = distFromFirstCorner[1] - (thirdCornerVec[1] - 2.0f * unskewF3D);
	distFromThirdCorner[2] = distFromFirstCorner[2] - (thirdCornerVec[2] - 2.0f * unskewF3D);
	distFromThirdCorner[3] = distFromFirstCorner[3] - (thirdCornerVec[3] - 2.0f * unskewF3D);

	distFromForthCorner[0] = distFromFirstCorner[0] - (forthCornerVec[0] - 3.0f * unskewF3D);
	distFromForthCorner[1] = distFromFirstCorner[1] - (forthCornerVec[1] - 3.0f * unskewF3D);
	distFromForthCorner[2] = distFromFirstCorner[2] - (forthCornerVec[2] - 3.0f * unskewF3D);
	distFromForthCorner[3] = distFromFirstCorner[3] - (forthCornerVec[3] - 3.0f * unskewF3D);

	distFromFifthCorner[0] = distFromFirstCorner[0] - (1.0f - 4.0f * unskewF3D);
	distFromFifthCorner[1] = distFromFirstCorner[1] - (1.0f - 4.0f * unskewF3D);
	distFromFifthCorner[2] = distFromFirstCorner[2] - (1.0f - 4.0f * unskewF3D);
	distFromFifthCorner[3] = distFromFirstCorner[3] - (1.0f - 4.0f * unskewF3D);

	/* Get Gradients */
	gradVec[0] = m_permTable[m_permTable[firstCornerVec[1] + m_permTable[firstCornerVec[2] + m_permTable[firstCornerVec[3]]]] + firstCornerVec[0]];
	gradVec[1] = m_permTable[m_permTable[firstCornerVec[1] + secondCornerVec[1] + m_permTable[firstCornerVec[2] + secondCornerVec[2] + m_permTable[firstCornerVec[3] + secondCornerVec[3]]]] + firstCornerVec[0] + secondCornerVec[0]];
	gradVec[2] = m_permTable[m_permTable[firstCornerVec[1] + thirdCornerVec[1] + m_permTable[firstCornerVec[2] + thirdCornerVec[2] + m_permTable[firstCornerVec[3] + thirdCornerVec[3]]]] + firstCornerVec[0] + thirdCornerVec[0]];
	gradVec[3] = m_permTable[m_permTable[firstCornerVec[1] + forthCornerVec[1] + m_permTable[firstCornerVec[2] + forthCornerVec[2] + m_permTable[firstCornerVec[3] + forthCornerVec[3]]]] + firstCornerVec[0] + forthCornerVec[0]];
	gradVec[4] = m_permTable[m_permTable[firstCornerVec[1] + 1 + m_permTable[firstCornerVec[2] + 1 + m_permTable[firstCornerVec[3] + 1]]] + firstCornerVec[0] + 1];

	/* Calculate noise contributions from each corner */
	noiseContrVec[0] = 0.5f - (distFromFirstCorner[0] * distFromFirstCorner[0]) - (distFromFirstCorner[1] * distFromFirstCorner[1]) - (distFromFirstCorner[2] * distFromFirstCorner[2]) - (distFromFirstCorner[3] * distFromFirstCorner[3]);
	noiseContrVec[1] = 0.5f - (distFromSecondCorner[0] * distFromSecondCorner[0]) - (distFromSecondCorner[1] * distFromSecondCorner[1]) - (distFromSecondCorner[2] * distFromSecondCorner[2]) - (distFromSecondCorner[3] * distFromSecondCorner[3]);
	noiseContrVec[2] = 0.5f - (distFromThirdCorner[0] * distFromThirdCorner[0]) - (distFromThirdCorner[1] * distFromThirdCorner[1]) - (distFromThirdCorner[2] * distFromThirdCorner[2]) - (distFromThirdCorner[3] * distFromThirdCorner[3]);
	noiseContrVec[3] = 0.5f - (distFromForthCorner[0] * distFromForthCorner[0]) - (distFromForthCorner[1] * distFromForthCorner[1]) - (distFromForthCorner[2] * distFromForthCorner[2]) - (distFromForthCorner[3] * distFromForthCorner[3]);
	noiseContrVec[4] = 0.5f - (distFromFifthCorner[0] * distFromFifthCorner[0]) - (distFromFifthCorner[1] * distFromFifthCorner[1]) - (distFromFifthCorner[2] * distFromFifthCorner[2]) - (distFromFifthCorner[3] * distFromFifthCorner[3]);

	/* This part could be optimized  */
	if (noiseContrVec[0] > 0) {
		noiseVec[0] = pow(noiseContrVec[0], 4.0f) * gradDot4D(gradVec[0], distFromFirstCorner[0], distFromFirstCorner[1], distFromFirstCorner[2], distFromFirstCorner[3]); //bContribution * bContribution * bContribution * bContribution * 
	}
	if (noiseContrVec[1] > 0) {
		noiseVec[1] = pow(noiseContrVec[1], 4.0f) * gradDot4D(gradVec[1], distFromSecondCorner[0], distFromSecondCorner[1], distFromSecondCorner[2], distFromSecondCorner[3]); //bContribution * bContribution * bContribution * bContribution * 
	}
	if (noiseContrVec[2] > 0) {
		noiseVec[2] = pow(noiseContrVec[2], 4.0f) * gradDot4D(gradVec[2], distFromThirdCorner[0], distFromThirdCorner[1], distFromThirdCorner[2], distFromThirdCorner[3]); //bContribution * bContribution * bContribution * bContribution * 
	}
	if (noiseContrVec[3] > 0) {
		noiseVec[3] = pow(noiseContrVec[3], 4.0f) * gradDot4D(gradVec[3], distFromForthCorner[0], distFromForthCorner[1], distFromForthCorner[2], distFromForthCorner[3]); //bContribution * bContribution * bContribution * bContribution * 
	}
	if (noiseContrVec[4] > 0) {
		noiseVec[3] = pow(noiseContrVec[4], 4.0f) * gradDot4D(gradVec[4], distFromFifthCorner[0], distFromFifthCorner[1], distFromFifthCorner[2], distFromFifthCorner[3]); //bContribution * bContribution * bContribution * bContribution * 
	}

	return 27.0f * (noiseVec[0] + noiseVec[1] + noiseVec[2] + noiseVec[3] + noiseVec[4]);
	//try (noiseVec[0] | noiseVec[1] | noiseVec[2] | noiseVec[3] | noiseVec[4]);
}

/* given hash, x, y values takes the dot product between pseudo-random gradient vector and a point around unit-circle (2d) */
float SimplexNoise::gradDot2D(const uint8_t hash, const float x, const float y) const
{
	const int newHash = hash & 7;
	return x * (float)m_grad2DTable[newHash][0] + y * (float)m_grad2DTable[newHash][1];
}

/* given hash, x, y, z values takes the dot product between pseudo-random gradient vector and a point around unit-cube (3d)
* For 3D Perlin reccomends to take 12 points rather than 16.
*/
float SimplexNoise::gradDot3D(const uint8_t hash, const float x, const float y, const float z) const
{
	const int newHash = hash & 15; //or 11, 
	return x * (float)m_grad3DTable[newHash][0] + y * (float)m_grad3DTable[newHash][1] + z * (float)m_grad3DTable[newHash][2];
}

/* given hash, x, y, z, t values takes  dot product between pseudo-random gradient vector and a point around unit-3-sphere (4d)
* Note: 4D usually used to generate animation on procedurally generated 3d content (like ocean, volumetric clouds, smoke simulations)
*/
float SimplexNoise::gradDot4D(const uint8_t hash, const float x, const float y, const float z, const float t) const
{
	const int newHash = hash & 31;
	return x * (float)m_grad4DTable[newHash][0] + y * (float)m_grad4DTable[newHash][1] + z * (float)m_grad4DTable[newHash][2] + t * (float)m_grad4DTable[newHash][3];
}

/* SETTERS */
void SimplexNoise::setOctaves(int octaves)
{
	m_octaves = octaves;
}
void SimplexNoise::setFrequency(float frequency)
{
	m_freq = frequency;
}
void SimplexNoise::setAmplitude(float amplitude)
{
	m_amp = amplitude;
}
void SimplexNoise::setSeed(int seed)
{
	m_seed = seed;
}
void SimplexNoise::setLacunarity(float lacunarity)
{
	m_lacunarity = lacunarity;
}
void SimplexNoise::setGain(float gain)
{
	m_gain = gain;
}


/* GETTERS */
int SimplexNoise::getOctaves() const
{
	return m_octaves;
}
float SimplexNoise::getFrequency() const
{
	return m_freq;
}
float SimplexNoise::getAmplitude() const
{
	return m_amp;
}
int SimplexNoise::getSeed() const
{
	return m_seed;
}
float SimplexNoise::getLacunarity() const
{
	return m_lacunarity;
}
float SimplexNoise::getGain() const
{
	return m_gain;
}


SimplexNoise::~SimplexNoise()
{
}
