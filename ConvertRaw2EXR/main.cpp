#include <ImfOutputFile.h>
#include <ImfChannelList.h>
#include <ImfArray.h>

#include <string>
#include <vector>
#include <sstream>
#include <fstream>

void encodeOpenExr(const std::string& fileName, double* image, int width, int height);

int main(int argc, char* argv[])
{
	if (argc < 5)
	{
		std::cout << "Please specify base file name, width, height, first and last frame!" << std::endl;
		return 0;
	}

	std::string baseFileName(argv[1]);

	int width = std::stoi(argv[2]);
	int height = std::stoi(argv[3]);

	int startFrame = std::stoi(argv[4]);
	int endFrame = std::stoi(argv[5]);

	std::vector<double> buffer(width * height * 3);

	for (int f = startFrame; f <= endFrame; ++f)
	{
		std::stringstream rawFileName;
		std::stringstream exrFileName;

		rawFileName << baseFileName << f << ".raw";
		exrFileName << baseFileName << f << ".exr";
		//1. read binary file
		std::ifstream file;

		file.open(rawFileName.str(), std::ios::binary);

		if (!file.is_open())
		{
			std::cout << "Error opening file: [ " << rawFileName.str() << " ]" << std::endl;
			continue;
		}

		file.read((char*)&buffer[0], width * height * 3 * (sizeof(double) / sizeof(char)));

		file.close();
		encodeOpenExr(exrFileName.str(), &buffer[0], width, height);

		rawFileName.clear();
		exrFileName.clear();
	}
}


//(almost straight from the documentation)
void encodeOpenExr(const std::string& fileName, double* image, int width, int height)
{
	Imath::Box2i displayWindow;
	Imath::Box2i dataWindow;

	Imf::Array2D<float> rPixels;
	Imf::Array2D<float> gPixels;
	Imf::Array2D<float> bPixels;

	rPixels.resizeErase(width, height);
	gPixels.resizeErase(width, height);
	bPixels.resizeErase(width, height);

	for (int row = 0; row < height; ++row)
	{
		for (int col = 0; col < width; ++col)
		{
			rPixels[row][col] = static_cast<float>(image[(row * width + col) * 3 + 0]);
			gPixels[row][col] = static_cast<float>(image[(row * width + col) * 3 + 1]);
			bPixels[row][col] = static_cast<float>(image[(row * width + col) * 3 + 2]);
		}
	}

	try
	{
		Imf::Header header(displayWindow, dataWindow);

		header.channels().insert("R", Imf::Channel(Imf::FLOAT));
		header.channels().insert("G", Imf::Channel(Imf::FLOAT));
		header.channels().insert("B", Imf::Channel(Imf::FLOAT));


		Imf::OutputFile file(fileName.c_str(), header);
		Imf::FrameBuffer frameBuffer;

		int dwWidth = dataWindow.max.x - dataWindow.min.x + 1;

		frameBuffer.insert("R",
			Imf::Slice(Imf::FLOAT,
			(char *)(&rPixels[0][0] - dataWindow.min.x - dataWindow.min.y * dwWidth),
			sizeof (rPixels[0][0]) * 1,
			sizeof (rPixels[0][0]) * dwWidth));

		frameBuffer.insert("G",
			Imf::Slice(Imf::FLOAT,
			(char *)(&gPixels[0][0] - dataWindow.min.x - dataWindow.min.y * dwWidth),
			sizeof (gPixels[0][0]) * 1,
			sizeof (gPixels[0][0]) * dwWidth));

		frameBuffer.insert("B",
			Imf::Slice(Imf::FLOAT,
			(char *)(&bPixels[0][0] - dataWindow.min.x - dataWindow.min.y * dwWidth),
			sizeof (bPixels[0][0]) * 1,
			sizeof (bPixels[0][0]) * dwWidth));

		file.setFrameBuffer(frameBuffer);
		file.writePixels(header.dataWindow().max.y - header.dataWindow().min.y + 1);
	}
	catch (Iex::BaseExc & e)
	{
		std::cerr << e.what() << std::endl;
	}
}