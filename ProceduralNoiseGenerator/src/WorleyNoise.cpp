#include "WorleyNoise.h"

/* CONSTANT VARIABLE INITIALIZATIONS */

uint8_t WorleyNoise::m_poissionTable[256] = { 4, 3, 1, 1, 1, 2, 4, 2, 2, 2, 5, 1, 0, 2, 1, 2, 2, 0, 4, 3, 2, 1, 2, 1, 3, 2, 2, 4, 2, 2, 5, 1, 2, 3,
2, 2, 2, 2, 2, 3, 2, 4, 2, 5, 3, 2, 2, 2, 5, 3, 3, 5, 2, 1, 3, 3, 4, 4, 2, 3, 0, 4, 2, 2, 2, 1, 3, 2,
2, 2, 3, 3, 3, 1, 2, 0, 2, 1, 1, 2, 2, 2, 2, 5, 3, 2, 3, 2, 3, 2, 2, 1, 0, 2, 1, 1, 2, 1, 2, 2, 1, 3,
4, 2, 2, 2, 5, 4, 2, 4, 2, 2, 5, 4, 3, 2, 2, 5, 4, 3, 3, 3, 5, 2, 2, 2, 2, 2, 3, 1, 1, 4, 2, 1, 3, 3,
4, 3, 2, 4, 3, 3, 3, 4, 5, 1, 4, 2, 4, 3, 1, 2, 3, 5, 3, 2, 1, 3, 1, 3, 3, 3, 2, 3, 1, 5, 5, 4, 2, 2,
4, 1, 3, 4, 1, 5, 3, 3, 5, 3, 4, 3, 2, 2, 1, 1, 1, 1, 1, 2, 4, 5, 4, 5, 4, 2, 1, 5, 1, 1, 2, 3, 3, 3,
2, 5, 2, 3, 3, 2, 0, 2, 1, 1, 4, 2, 1, 3, 2, 1, 2, 2, 3, 2, 5, 5, 3, 4, 5, 5, 2, 4, 4, 5, 3, 2, 2, 2,
1, 4, 2, 3, 3, 4, 2, 5, 4, 2, 4, 2, 2, 2, 4, 5, 3, 2 };

const float WorleyNoise::DENSITY_ADJUSTMENT = 0.398150f; //this may not be correct for every single distance funct. 
const float WorleyNoise::INV_DENSITY_ADJ = (1.0f / DENSITY_ADJUSTMENT);
const float WorleyNoise::MINKOWSKI_INV = (1.0f / MINKOWSKI_CONST);
const float WorleyNoise::SCALE_CONSTANT = 0.000000000232831f; //(1.0f / 4294967296.0f)



WorleyNoise::WorleyNoise(int octaves, float frequency, float amplitude, int seed, float lacunarity, float gain, distanceFunctions funct, noiseTypes noiseType )
{
	m_seed = seed;
	m_octaves = octaves;
	m_freq = frequency;
	m_amp = amplitude;
	m_lacunarity = lacunarity;
	m_gain = gain;
	m_func = funct;
	m_noiseType = noiseType;
	pickDistanceFunct();
	srand(m_seed);

	std::random_shuffle(std::begin(m_poissionTable), std::end(m_poissionTable)); //shuffle the distribution table
}

/* The most naive way of checking closes neighbors is to test against every boundary shape (cubes in 3d, squares in 2d)
* But that's way too slow, Worley introduced a more efficient model (http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.95.412&rep=rep1&type=pdf)
* Instead, we can just check if it's possible for neighbor shapes to be contributors by examining the combinations of the sum
of the squared distances from the shapes lower and upper corners
*/

/* fractal brownian motion 2D */
float WorleyNoise::worleyNoise2D(float x, float y)
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

/* fractal brownian motion 3D */
float WorleyNoise::worleyNoise3D(float x, float y, float z)
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
/* No smoothing or anti-aliasing done on the noise, just takes the abs of the noise */
float WorleyNoise::turbulanceWorleyNoise2D(float x, float y)
{
	float tempAmp = m_amp;
	float tempFreq = m_freq;
	float pixNoise = 0.0f;

	for (int i = 0; i < m_octaves; ++i) {
		float noise = noise2D(x * tempFreq, y * tempFreq);
		pixNoise += fabs(noise) * tempAmp;
		tempFreq *= m_lacunarity;
		tempAmp *= m_gain;
	}
	return pixNoise;
}

/* No smoothing or anti-aliasing done on the noise, just takes the abs of the noise */
float WorleyNoise::turbulanceWorleyNoise3D(float x, float y, float z)
{
	float tempAmp = m_amp;
	float tempFreq = m_freq;
	float pixNoise = 0.0f;

	for (int i = 0; i < m_octaves; ++i) {
		float noise = noise3D(x * tempFreq, y * tempFreq, z * tempFreq);
		pixNoise += fabs(noise) * tempAmp;
		tempFreq *= m_lacunarity;
		tempAmp *= m_gain;
	}
	return pixNoise;
}


/* The most naive way of checking closes neighbors is to test against every boundary shape (cubes in 3d, squares in 2d)
* But that's way too slow, Worley introduced a more efficient model (http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.95.412&rep=rep1&type=pdf)
* Instead, we can just check if it's possible for neighbor shapes to be contributors by examining the combinations of the sum
of the squared distances from the shapes lower and upper corners
*/

float WorleyNoise::noise2D(float x, float y)
{
	/*upper & lower corner vals for x and y. [0] -> x, [1] -> y */
	float upperCorners[2], lowerCorners[2]; //float xUpperCorner, yUpperCorner, xLowerCorner, yLowerCorner;
	/* Reset feature points to the highest possible float val */
	std::fill_n(m_featurePts, NUM_OF_FEATURES, FLT_MAX);

	/* Apply density adjustment on coordinates */
	float adjustedCords[2] = { DENSITY_ADJUSTMENT * x, DENSITY_ADJUSTMENT * y };
	int intCords[2] = { NoiseHelperFunctions::fastFloor(adjustedCords[0]), NoiseHelperFunctions::fastFloor(adjustedCords[1]) };

	/* Add the first adjusted coordinates to m_delta2D */
	addSamples2D(intCords[0], intCords[1], adjustedCords);

	/* Check if it's possible for neighbor shape to be a contributor. */
	lowerCorners[0] = adjustedCords[0] - intCords[0];
	lowerCorners[1] = adjustedCords[1] - intCords[1];
	upperCorners[0] = (1.0f - lowerCorners[0]) * (1.0f - lowerCorners[0]);
	upperCorners[1] = (1.0f - lowerCorners[1]) * (1.0f - lowerCorners[1]);
	lowerCorners[0] *= lowerCorners[0]; //x - 1 + xLowerCorner
	lowerCorners[1] *= lowerCorners[1];

	/* Test if any of the 4 Vertices would be in feature coordinates.  */
	if (lowerCorners[0] < m_featurePts[NUM_OF_FEATURES - 1])	addSamples2D(intCords[0] - 1, intCords[1], adjustedCords);
	if (lowerCorners[1] < m_featurePts[NUM_OF_FEATURES - 1])	addSamples2D(intCords[0], intCords[1] - 1, adjustedCords);
	if (upperCorners[0] < m_featurePts[NUM_OF_FEATURES - 1])	addSamples2D(intCords[0] + 1, intCords[1], adjustedCords);
	if (upperCorners[1] < m_featurePts[NUM_OF_FEATURES - 1])	addSamples2D(intCords[0], intCords[1] + 1, adjustedCords);

	/* Test if any of the 4 edges would be in feature coordinates. */
	if (lowerCorners[0] + lowerCorners[1] < m_featurePts[NUM_OF_FEATURES - 1])	addSamples2D(intCords[0] - 1, intCords[1] - 1, adjustedCords);
	if (lowerCorners[0] + upperCorners[1] < m_featurePts[NUM_OF_FEATURES - 1])	addSamples2D(intCords[0] - 1, intCords[1] + 1, adjustedCords);
	if (upperCorners[0] + upperCorners[1] < m_featurePts[NUM_OF_FEATURES - 1])	addSamples2D(intCords[0] + 1, intCords[1] + 1, adjustedCords);
	if (upperCorners[0] + lowerCorners[1] < m_featurePts[NUM_OF_FEATURES - 1])	addSamples2D(intCords[0] + 1, intCords[1] - 1, adjustedCords);

	//Scale values back to their original scale
	for (int i = 0; i < NUM_OF_FEATURES; ++i) {
		m_featurePts[i] = NoiseHelperFunctions::fastSqrt(m_featurePts[i]) * INV_DENSITY_ADJ; //use fastsqrt later
		m_delta2D[i][0] *= INV_DENSITY_ADJ;
		m_delta2D[i][1] *= INV_DENSITY_ADJ;
	}


	//float fsum = m_featurePts[0] + m_featurePts[1];
	//float fdelta = m_featurePts[0] - m_featurePts[1];
	return noiseType();
}

/* 3D requires to check 8 vertices, 12 edges, and 6 faces. */
float WorleyNoise::noise3D(float x, float y, float z)
{
	float upperCorners[3], lowerCorners[3]; //float xUpperCorner, yUpperCorner, zUpperCorner, xLowerCorner, yLowerCorner, zLowerCorner;

	/* Reset feature points to the highest possible float val */
	std::fill_n(m_featurePts, NUM_OF_FEATURES, FLT_MAX);

	/* Apply density adjustment on coordinates */
	float adjustedCords[3] = { DENSITY_ADJUSTMENT * x, DENSITY_ADJUSTMENT * y, DENSITY_ADJUSTMENT * z };
	int intCords[3] = { NoiseHelperFunctions::fastFloor(adjustedCords[0]), NoiseHelperFunctions::fastFloor(adjustedCords[1]), NoiseHelperFunctions::fastFloor(adjustedCords[2]) };

	/* Add the first adjusted coordinates to m_delta3D */
	addSamples3D(intCords[0], intCords[1], intCords[2], adjustedCords);

	/* Get upper and lower corners of the cube to check if it's possible for neighbor shape to be a contributor. */
	lowerCorners[0] = adjustedCords[0] - intCords[0];
	lowerCorners[1] = adjustedCords[1] - intCords[1];
	lowerCorners[2] = adjustedCords[2] - intCords[2];

	upperCorners[0] = (1.0f - lowerCorners[0]) * (1.0f - lowerCorners[0]);
	upperCorners[1] = (1.0f - lowerCorners[1]) * (1.0f - lowerCorners[1]);
	upperCorners[2] = (1.0f - lowerCorners[2]) * (1.0f - lowerCorners[2]);

	lowerCorners[0] *= lowerCorners[0];
	lowerCorners[1] *= lowerCorners[1];
	lowerCorners[2] *= lowerCorners[2];

	/* Test if any of the 6 faces would be in feature coordinates.  */
	if (lowerCorners[0] < m_featurePts[NUM_OF_FEATURES - 1])	addSamples3D(intCords[0] - 1, intCords[1], intCords[2], adjustedCords);
	if (lowerCorners[1] < m_featurePts[NUM_OF_FEATURES - 1])	addSamples3D(intCords[0], intCords[1] - 1, intCords[2], adjustedCords);
	if (lowerCorners[2] < m_featurePts[NUM_OF_FEATURES - 1])	addSamples3D(intCords[0], intCords[1], intCords[2] - 1, adjustedCords);
	if (upperCorners[0] < m_featurePts[NUM_OF_FEATURES - 1])	addSamples3D(intCords[0] + 1, intCords[1], intCords[2], adjustedCords);
	if (upperCorners[1] < m_featurePts[NUM_OF_FEATURES - 1])	addSamples3D(intCords[0], intCords[1] + 1, intCords[2], adjustedCords);
	if (upperCorners[2] < m_featurePts[NUM_OF_FEATURES - 1])	addSamples3D(intCords[0], intCords[1], intCords[2] + 1, adjustedCords);


	/* Test if any of the 12 edges would be in feature coordinates. */
	if (lowerCorners[0] + lowerCorners[1] < m_featurePts[NUM_OF_FEATURES - 1])	addSamples3D(intCords[0] - 1, intCords[1] - 1, intCords[2], adjustedCords);
	if (lowerCorners[0] + lowerCorners[2] < m_featurePts[NUM_OF_FEATURES - 1])	addSamples3D(intCords[0] - 1, intCords[1], intCords[2] - 1, adjustedCords);
	if (lowerCorners[1] + lowerCorners[2] < m_featurePts[NUM_OF_FEATURES - 1])	addSamples3D(intCords[0], intCords[1] - 1, intCords[2] - 1, adjustedCords);

	if (upperCorners[0] + upperCorners[1] < m_featurePts[NUM_OF_FEATURES - 1])	addSamples3D(intCords[0] + 1, intCords[1] + 1, intCords[2], adjustedCords);
	if (upperCorners[0] + upperCorners[2] < m_featurePts[NUM_OF_FEATURES - 1])	addSamples3D(intCords[0] + 1, intCords[1], intCords[2] + 1, adjustedCords);
	if (upperCorners[1] + upperCorners[2] < m_featurePts[NUM_OF_FEATURES - 1])	addSamples3D(intCords[0], intCords[1] + 1, intCords[2] + 1, adjustedCords);

	if (lowerCorners[0] + upperCorners[1] < m_featurePts[NUM_OF_FEATURES - 1])	addSamples3D(intCords[0] - 1, intCords[1] + 1, intCords[2], adjustedCords);
	if (lowerCorners[0] + upperCorners[2] < m_featurePts[NUM_OF_FEATURES - 1])	addSamples3D(intCords[0] - 1, intCords[1], intCords[2] + 1, adjustedCords);
	if (lowerCorners[1] + upperCorners[2] < m_featurePts[NUM_OF_FEATURES - 1])	addSamples3D(intCords[0], intCords[1] - 1, intCords[2] + 1, adjustedCords);

	if (upperCorners[0] + lowerCorners[1] < m_featurePts[NUM_OF_FEATURES - 1])	addSamples3D(intCords[0] + 1, intCords[1] - 1, intCords[2], adjustedCords);
	if (upperCorners[0] + lowerCorners[2] < m_featurePts[NUM_OF_FEATURES - 1])	addSamples3D(intCords[0] + 1, intCords[1], intCords[2] - 1, adjustedCords);
	if (upperCorners[1] + lowerCorners[2] < m_featurePts[NUM_OF_FEATURES - 1])	addSamples3D(intCords[0], intCords[1] + 1, intCords[2] - 1, adjustedCords);

	/* Test if any of the 8 vertices would be in feature coordinates. */
	if (lowerCorners[0] + lowerCorners[1] + lowerCorners[2] < m_featurePts[NUM_OF_FEATURES - 1])	addSamples3D(intCords[0] - 1, intCords[1] - 1, intCords[2] - 1, adjustedCords);
	if (lowerCorners[0] + lowerCorners[1] + upperCorners[2] < m_featurePts[NUM_OF_FEATURES - 1])	addSamples3D(intCords[0] - 1, intCords[1] - 1, intCords[2] + 1, adjustedCords);
	if (lowerCorners[0] + upperCorners[1] + lowerCorners[2] < m_featurePts[NUM_OF_FEATURES - 1])	addSamples3D(intCords[0] - 1, intCords[1] + 1, intCords[2] - 1, adjustedCords);
	if (lowerCorners[0] + upperCorners[1] + upperCorners[2] < m_featurePts[NUM_OF_FEATURES - 1])	addSamples3D(intCords[0] - 1, intCords[1] + 1, intCords[2] + 1, adjustedCords);

	if (upperCorners[0] + lowerCorners[1] + lowerCorners[2] < m_featurePts[NUM_OF_FEATURES - 1])	addSamples3D(intCords[0] + 1, intCords[1] - 1, intCords[2] - 1, adjustedCords);
	if (upperCorners[0] + lowerCorners[1] + upperCorners[2] < m_featurePts[NUM_OF_FEATURES - 1])	addSamples3D(intCords[0] + 1, intCords[1] - 1, intCords[2] + 1, adjustedCords);
	if (upperCorners[0] + upperCorners[1] + lowerCorners[2] < m_featurePts[NUM_OF_FEATURES - 1])	addSamples3D(intCords[0] + 1, intCords[1] + 1, intCords[2] - 1, adjustedCords);
	if (upperCorners[0] + upperCorners[1] + upperCorners[2] < m_featurePts[NUM_OF_FEATURES - 1])	addSamples3D(intCords[0] + 1, intCords[1] + 1, intCords[2] + 1, adjustedCords);

	//Scale values back to their original scale
	for (int i = 0; i < NUM_OF_FEATURES; ++i) {
		m_featurePts[i] = fastSqrt(m_featurePts[i]) * INV_DENSITY_ADJ;
		m_delta3D[i][0] *= INV_DENSITY_ADJ;
		m_delta3D[i][1] *= INV_DENSITY_ADJ;
		m_delta3D[i][2] *= INV_DENSITY_ADJ;
	}


	/* This part should be a mode to select in real time. */
	//float fsum = m_featurePts[0] + m_featurePts[1];
	//float fdelta = m_featurePts[1] - m_featurePts[0];

	return noiseType();
}

/* Sort square of samples into the current
* We'll be dealing with the int parts of coordinates.
*/
void  WorleyNoise::addSamples2D(int x, int y, const float(&newPos)[2])
{
	float featurePoints[2];
	float deltas[2];
	float distance;
	uint8_t fCount;
	unsigned long seed = 702395077 * x + 915488749 * y;
	/* Get the size of feature coordinates in this square. */

	//cout << (seed >> 24) << endl;
	fCount = m_poissionTable[seed >> 24]; // no need for 256
	NoiseHelperFunctions::LCG(seed); //churn them seeds
	for (int i = 0; i < fCount; ++i) {
		//if you want to store ids, store seed before you churn it.
		NoiseHelperFunctions::LCG(seed);

		/* get Feature points for x and y and scale it to be between 0 and 1 */
		featurePoints[0] = (seed + 0.5f) * SCALE_CONSTANT;
		//cout << xFeaturePt << endl;
		NoiseHelperFunctions::LCG(seed);
		featurePoints[1] = (seed + 0.5f) * SCALE_CONSTANT;
		NoiseHelperFunctions::LCG(seed);

		/* Calculate deltas of feature points */
		deltas[0] = x + featurePoints[0] - newPos[0];
		deltas[1] = y + featurePoints[1] - newPos[1];

		distance = (this->*m_distanceFunct2D)(deltas);

		/* If distance higher than the largest order of feature point then it's not in the circle, don't bother inserting it to delta2D */
		if (distance < m_featurePts[NUM_OF_FEATURES - 1]) {
			int index = NUM_OF_FEATURES;
			while (index > 0 && (distance < m_featurePts[index - 1])) { index--; };

			for (int i = NUM_OF_FEATURES - 2; i >= index; --i)
			{
				m_featurePts[i + 1] = m_featurePts[i];
				m_delta2D[i + 1][0] = m_delta2D[i][0];
				m_delta2D[i + 1][1] = m_delta2D[i][1];
			}
			// Insert the new point's information into the list.
			m_featurePts[index] = distance;
			m_delta2D[index][0] = deltas[0];
			m_delta2D[index][1] = deltas[1];
		}
	}
}

void  WorleyNoise::addSamples3D(int x, int y, int z, const float(&newPos)[3])
{
	float featurePoints[3];
	float deltas[3];
	float distance;
	uint8_t fCount;

	unsigned long seed = 702395077 * x + 915488749 * y + 2120969693 * z;
	/* Get the size of feature coordinates in this square. */
	fCount = m_poissionTable[seed >> 24];
	LCG(seed); //churn them seeds
	for (int i = 0; i < fCount; ++i) {
		//if you want to store ids, store seed before you churn it.
		LCG(seed);

		/* get Feature points for x and y and scale it to be between 0 and 1 */
		featurePoints[0] = (seed + 0.5f) * SCALE_CONSTANT;
		LCG(seed);
		featurePoints[1] = (seed + 0.5f) * SCALE_CONSTANT;
		LCG(seed);
		featurePoints[2] = (seed + 0.5f) * SCALE_CONSTANT;
		LCG(seed);

		/* Calculate deltas of feature points */
		deltas[0] = x + featurePoints[0] - newPos[0];
		deltas[1] = y + featurePoints[1] - newPos[1];
		deltas[2] = z + featurePoints[2] - newPos[2];

		distance = (this->*m_distanceFunct3D)(deltas);

		/* If distance higher than the largest order of feature point then it's not in the circle, don't bother inserting it to delta2D */
		if (distance < m_featurePts[NUM_OF_FEATURES - 1]) {
			int index = NUM_OF_FEATURES;
			while (index > 0 && (distance < m_featurePts[index - 1])) { index--; };
			for (int i = NUM_OF_FEATURES - 2; i >= index; --i) {
				m_featurePts[i + 1] = m_featurePts[i];
				m_delta3D[i + 1][0] = m_delta3D[i][0];
				m_delta3D[i + 1][1] = m_delta3D[i][1];
				m_delta3D[i + 1][2] = m_delta3D[i][2];
			}
			// Insert the new point's information into the list.
			m_featurePts[index] = distance;
			m_delta3D[index][0] = deltas[0];
			m_delta3D[index][1] = deltas[1];
			m_delta3D[index][2] = deltas[2];
		}
	}
}

/* Manipulating Feature points could yield wildly different results.
* These are only few I implemented, experimenting with new types higly encouraged!
*/
float WorleyNoise::noiseType()
{
	switch (NUM_OF_FEATURES) {
	case 1:
		return m_featurePts[0];
	case 2:
		switch (m_noiseType) {
		case FIRST_ORDER:
			return m_featurePts[0];
		case SECOND_ORDER:
			return m_featurePts[1];
		case SUM_F1F2:
			return m_featurePts[0] + m_featurePts[1];
		case  DELTA_F2F1:
			return m_featurePts[1] - m_featurePts[0];
		case MULTIPLY_F1F2:
			return m_featurePts[0] * m_featurePts[1];
		case DIVIDE_F2F1:
			return m_featurePts[1] / m_featurePts[0];
		default:
			throw std::invalid_argument("Error: Declared noise type is not supported for current number of Feature points!");
		}
	default: //or case 3 if you have case 4:
		switch (m_noiseType) {
		case FIRST_ORDER:
			return m_featurePts[0];
		case SECOND_ORDER:
			return m_featurePts[1];
		case THIRD_ORDER:
			return m_featurePts[2];

		case SUM_F1F2:
			return m_featurePts[0] + m_featurePts[1];
		case SUM_F1F3:
			return m_featurePts[0] + m_featurePts[2];
		case SUM_F2F3:
			return m_featurePts[1] + m_featurePts[2];
		case SUM_F1F2F3:
			return m_featurePts[0] + m_featurePts[1] + m_featurePts[2];

		case DELTA_F2F1:
			return m_featurePts[1] - m_featurePts[0];
		case DELTA_F3F2:
			return  m_featurePts[2] - m_featurePts[1];
		case DELTA_F3F1:
			return  m_featurePts[2] - m_featurePts[0];
		case DELTA_F3F2F1:
			return  m_featurePts[2] - m_featurePts[1] - m_featurePts[0];

		case MULTIPLY_F1F2:
			return m_featurePts[0] * m_featurePts[1];
		case MULTIPLY_F1F3:
			return m_featurePts[0] * m_featurePts[2];
		case MULTIPLY_F2F3:
			return m_featurePts[1] * m_featurePts[2];
		case MULTIPLY_F1F2F3:
			return  m_featurePts[2] * m_featurePts[1] * m_featurePts[0];	// DIVIDE_F2F1, DIVIDE_F3F2, DIVIDE_F3F1, DIVIDE_F3F2F1

		case DIVIDE_F2F1:
			return m_featurePts[1] / m_featurePts[0];
		case DIVIDE_F3F1:
			return m_featurePts[2] / m_featurePts[0];
		case DIVIDE_F3F2:
			return m_featurePts[2] / m_featurePts[1];
		case DIVIDE_F3F2F1:
			return m_featurePts[2] / m_featurePts[1] / m_featurePts[0];
		default:
			throw std::invalid_argument("Error: Declared noise type is not supported for current number of Feature points!");
		}
	}

}

/* Originally I wanted to use a function pointer that points to the distance function specified by m_func, I was bit cleaner but much more slower than using switches */
void WorleyNoise::pickDistanceFunct()
{
	switch (m_func) {
	case EUCLIDEAN_SQUARED:
		m_distanceFunct2D = &WorleyNoise::euclideanDistanceSquared;
		m_distanceFunct3D = &WorleyNoise::euclideanDistanceSquared;
		break;

	case QUADRATIC:
		m_distanceFunct2D = &WorleyNoise::quadraticDistance;
		m_distanceFunct3D = &WorleyNoise::quadraticDistance;
		break;

	case EUCLIDEAN:
		m_distanceFunct2D = &WorleyNoise::euclideanDistance;
		m_distanceFunct3D = &WorleyNoise::euclideanDistance;
		break;

	case MANHATTAN:
		m_distanceFunct2D = &WorleyNoise::manhattanDistance;
		m_distanceFunct3D = &WorleyNoise::manhattanDistance;
		break;

	case MINKOWSKI:
		m_distanceFunct2D = &WorleyNoise::minkowskiDistance;
		m_distanceFunct3D = &WorleyNoise::minkowskiDistance;
		break;

	case CITYBLOCK:
		m_distanceFunct2D = &WorleyNoise::cityBlockDistance;
		m_distanceFunct3D = &WorleyNoise::cityBlockDistance;
		break;
	default:
		throw std::invalid_argument("Error: Undefined Distance Function!");
	}
}

/* SETTERS */
void WorleyNoise::setOctaves(int octaves)
{
	m_octaves = octaves;
}
void WorleyNoise::setFrequency(float frequency)
{
	m_freq = frequency;
}
void WorleyNoise::setAmplitude(float amplitude)
{
	m_amp = amplitude;
}
void WorleyNoise::setSeed(int seed)
{
	m_seed = seed;
}
void WorleyNoise::setLacunarity(float lacunarity)
{
	m_lacunarity = lacunarity;
}
void WorleyNoise::setGain(float gain)
{
	m_gain = gain;
}
void WorleyNoise::setDistanceFunction(distanceFunctions distFunct)
{
	m_func = distFunct;
}
void WorleyNoise::setNoiseType(noiseTypes noiseType)
{
	m_noiseType = noiseType;
}

/* GETTERS */
int WorleyNoise::getOctaves() const
{
	return m_octaves;
}
float WorleyNoise::getFrequency() const
{
	return m_freq;
}
float WorleyNoise::getAmplitude() const
{
	return m_amp;
}
int WorleyNoise::getSeed() const
{
	return m_seed;
}
float WorleyNoise::getLacunarity() const
{
	return m_lacunarity;
}
float WorleyNoise::getGain() const
{
	return m_gain;
}
WorleyNoise::distanceFunctions WorleyNoise::getDistanceFunction() const
{
	return m_func;
}
WorleyNoise::noiseTypes WorleyNoise::getNoiseType() const
{
	return m_noiseType;
}

WorleyNoise::~WorleyNoise()
{
}
