#pragma once
#pragma warning( disable : 4996 )
#include <iostream>
#include <exception>
#include "MapGenerator.h"



using namespace std;



/* 
 Here's a short info on what each attribute does when generating a noise:

* Octaves: Number of interpolation pass done on the noise, it effectively controls level of detail. Higher the octave, more detailed the image looks.

* Frequency: Determines number of terrain features and size of these features. Decreasing it will decrease number of features and increase size of the features.
There might not be enough pixels to show the details, in fact if frequency is > resolution of the heightmap extra details won't have any affect.

* Amplitude: Amplitude determines the max value for the noise on 3d, amplitude can affect the height of the terrain.
For each octave, amplitude shrinks, when amplitudes are too small, there is no point adding more noise functions since it won't have any significant effect.

* Seed: random seed to pass to srand() which will effect the shuffle pattern (same seed will result in same shuffle pattern), you might want to play around with this value to find interesting patterns.

* Lacunarity: Determines how frequency grows. Usually for each octave, frequency doubles in value.

* Gain: Determines how much the amplitude shrink. Usually for each octave, amplitude halved in value.

* You can see the results of how gain and lacunarity affect signals here: http://freespace.virgin.net/hugo.elias/models/m_perlin.htm

To change the color set for the map, please check the MapGenerator class, simple overview of what function/class does located in .h, nitty-gritty details located in .cpps
*/


/* Main is only there to show how each function works and called, it's entirely possible to generate multiple noises with the same generator, instead of selecting one at a time.
* Noise Classes are loosely-coupled on purpose so that it's easy to import them to your project, Note that ATM all the noises are generated on CPU, if you want to use them in real-time applications, you're better off  implementing them as shaders on GPU instead.
*/

/* TODO: 
 *	write the --help 
 *  Implement a simple cross platform GUI to get rid of this clunky text based input acquisition.
 */



int main()
{

	MapGenerator::Noise noiseType;
	WorleyNoise::distanceFunctions distanceFuncType;

	char noiseTypeChar;
	char distanceFuncTypeChar;
	char defaultValueInput;
	bool areValuesDefault;

	int hgrid, vgrid;
	float frequency;
	float octaves;
	float amplitude;
	int seed;
	float lacunarity;
	float gain;



	cout << "Enter Map dimensions: (integer width and height seperated by space) ";
	try {
		cin >> hgrid;
		cin >> vgrid;
	}
	catch (logic_error) {
		cout << "Error! You have entered an invalid input!" << endl;
	}
	
	MapGenerator generator(hgrid, vgrid, "OUTPUT");


	cout << "Please Enter Noise Type: (P for Perlin Noise, S for Simplex Noise, W for Worley Noise(Cell Noise) ";
	cin >> noiseTypeChar;

	//16 112 p  | 19 115 s | 23 119 w
	if (noiseTypeChar == 80 || noiseTypeChar == 112) {
		noiseType = MapGenerator::Noise::PERLIN;
	}
	else if (noiseTypeChar == 83 || noiseTypeChar == 115) {
		noiseType = MapGenerator::Noise::SIMPLEX;

	}
	else if (noiseTypeChar == 87 || noiseTypeChar == 119) {
		noiseType = MapGenerator::Noise::WORLEY;
	}
	else {
		cout << "Error! Invalid Noise Type please refer to --help to see supported noise types" << endl;
		exit(-1);
	}


	cout << "Would you like to use default octaves, frequency, amplitude, seed, lacunarity, gain values? (Y/N): ";
	cin >> defaultValueInput;

	if (defaultValueInput == 89 || defaultValueInput == 121) {
		areValuesDefault = true;
	}
	else if (defaultValueInput == 78 || defaultValueInput == 110) {
		areValuesDefault = false;
	}
	else {
		cout << "Error! Expected Y or N" << endl;
		exit(-1);
	}


	if (!areValuesDefault) {
		cout << "Please Enter octaves, frequency, amplitude, seed, lacunarity, gain (seperated by space): ";
		try {
			cin >> octaves;
			cin >> frequency;
			cin >> amplitude;
			cin >> seed;
			cin >> lacunarity;
			cin >> gain;
		}
		catch (logic_error) {
			cout << "Error! You have entered an invalid input!" << endl;
			exit(-1);
		}
	}
	else {
		frequency = 1.0f / (float)hgrid * 5;
		octaves = 17;
		amplitude = 3.5f;
		seed = 7;
		lacunarity = 2.0f;
		gain = 0.60f;

	}


	//Ask which distance function to be used if Worley Noise is selected
	if (noiseType == MapGenerator::Noise::WORLEY) {
		cout << "Please Enter which Distance Function to use for Worley: (Q for Quadratic, E for Euclidean, S for Euclidean Square, M for Manhattan, W for Minkowski, C for Cityblock) ";
		cin >> distanceFuncTypeChar;
		if (distanceFuncTypeChar == 81 || distanceFuncTypeChar == 113) {
			distanceFuncType = WorleyNoise::distanceFunctions::QUADRATIC;
		}
		else if (distanceFuncTypeChar == 69 || distanceFuncTypeChar == 101) {
			distanceFuncType = WorleyNoise::distanceFunctions::EUCLIDEAN;
		}
		else if (distanceFuncTypeChar == 83 || distanceFuncTypeChar == 115) {
			distanceFuncType = WorleyNoise::distanceFunctions::EUCLIDEAN_SQUARED;
		}
		else if (distanceFuncTypeChar == 77 || distanceFuncTypeChar == 109) {
			distanceFuncType = WorleyNoise::distanceFunctions::MANHATTAN;
		}
		else if (distanceFuncTypeChar == 87 || distanceFuncTypeChar == 119) {
			distanceFuncType = WorleyNoise::distanceFunctions::MINKOWSKI;
		}
		else if (distanceFuncTypeChar == 67 || distanceFuncTypeChar == 99) {
			distanceFuncType = WorleyNoise::distanceFunctions::MINKOWSKI;
		}
		else {
			cout << "Error! Invalid Distance Function! please refer to --help to see supported distance functions" << endl;
			exit(-1);
		}
	}
	

	//Noise funct created
	if (noiseType == MapGenerator::Noise::PERLIN) {
		generator.createPerlinNoise(octaves, frequency, amplitude, seed, lacunarity, gain);


	}
	else if (noiseType == MapGenerator::Noise::SIMPLEX) {
		generator.createSimplexNoise(octaves, frequency, amplitude, seed, lacunarity, gain);

	}
	else if (noiseType == MapGenerator::Noise::WORLEY) {
		generator.createWorleyNoise(octaves, frequency, amplitude, seed, lacunarity, gain, distanceFuncType);
	}


	generator.fillMapWithNoise2D(noiseType);
	generator.printGrayscaleImage("MAPTEST_GRAYSCALE.bmp");
	generator.printColoredImage("MAPTEST_COLORED.bmp");
	generator.flushMap();

	cout << "Textures are Generated!" << endl;
	return 0;

	/*
	hgrid = vgrid = 500;

	float frequency = 1.0f / (float)hgrid * 5;
	float octaves = 17;
	float amplitude = 3.5f;

	MapGenerator generator(hgrid, vgrid, "MAPTEST");
	generator.createPerlinNoise(octaves, frequency, amplitude, 666, 2.0f, 0.60f);
	generator.createSimplexNoise(octaves, frequency, amplitude, 666, 2.0f, 0.65f);
	generator.createWorleyNoise(octaves, frequency, amplitude, 666, 2.0f, 0.65f);
	generator.printNoiseInfo();

	generator.fillMapWithNoise2D(MapGenerator::Noise::PERLIN);
	generator.printGrayscaleImage("MAPTEST.bmp");
	generator.printColoredImage("MAPTEST_COLORED.bmp");
	generator.flushMap();

	generator.fillMapWithNoise2D(MapGenerator::Noise::SIMPLEX);
	generator.printGrayscaleImage("MAPTEST_SIMPLEX.bmp");
	generator.printColoredImage("MAPTEST_COLORED_SIMPLEX.bmp");
	generator.flushMap();

	generator.fillMapWithNoise2D(MapGenerator::Noise::WORLEY);
	generator.printGrayscaleImage("MAPTEST_WORLEY.bmp");
	generator.printColoredImage("MAPTEST_COLORED_WORLEY.bmp");
	generator.flushMap();

	cout << "Finished!" << endl;
	system("pause");
	return 0;
	*/
}

