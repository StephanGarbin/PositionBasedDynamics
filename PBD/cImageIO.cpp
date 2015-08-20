#include "cImageIO.h"
#include "lodepng.h"
#include <iostream>
#include <sstream>
#include <fstream>

void saveAsRaw(double* data, int width, int height, int frame, const std::string& fileName);
void saveAsPng(double* data, int width, int height, int frame, const std::string& fileName);

void saveImage(int width, int height, int frame, const std::string& fileName)
{
	double *image;

	/* Allocate our buffer for the image */
	if ((image = (double*)malloc(3 * width*height*sizeof(double))) == NULL) {
		fprintf(stderr, "Failed to allocate memory for image\n");
	}

	/* Copy the image into our buffer */

	GLenum err = glGetError();
	if (err != GL_NO_ERROR)
	{
		std::cerr << "GL Error (Before): " << err << std::endl;
	}

	glReadBuffer(GL_FRONT);
	err = glGetError();
	if (err != GL_NO_ERROR)
	{
		std::cerr << "GL Error (Read Buffer): " << err << std::endl;
	}
	glReadPixels(0, 0, width, height, GL_RGB, GL_DOUBLE, image);

	err = glGetError();
	if (err != GL_NO_ERROR)
	{
		std::cerr << "GL Error (Read Pixels: " << err << std::endl;
	}

	saveAsPng(image, width, height, frame, fileName);

	free(image);
}

void saveAsRaw(double* data, int width, int height, int frame, const std::string& fileName)
{
	std::stringstream completeFileName;
	completeFileName << "images/" << fileName << "_" << frame << ".txt";

	std::ofstream file;

	file.open(completeFileName.str().c_str(), std::ios::binary);

	if (!file.is_open())
	{
		std::cout << "Error: Could not open file: [ " << completeFileName.str() << " ] " << std::endl;
	}

	for (int i = 0; i < width * height * 3; ++i)
	{
		file << data[i];
	}

	file.close();
	file.clear();
	std::cout << "Written raw file: " << completeFileName.str() << std::endl;

}

void saveAsPng(double* data, int width, int height, int frame, const std::string& fileName)
{
	std::stringstream completeFileName;
	completeFileName << "images/" << fileName << "_" << frame << ".png";

	std::vector<unsigned char> image(width * height * 3);
	for (int i = 0; i < image.size(); ++i)
	{
		if (data[i] != 0.0)
		{
			std::cout << data[i] << ",";
		}
		image[i] = static_cast<unsigned char>(data[i] * 255);
	}

	//Encode the image
	unsigned error = lodepng::encode(completeFileName.str(), &image[0], width, height, LodePNGColorType::LCT_RGB);

	//if there's an error, display it
	if (error) std::cout << "encoder error " << error << ": " << lodepng_error_text(error) << std::endl;
}
