#include "AbcReader.h"

#include <vector>
#include <string>

int main(int argc, char* argv[])
{
	std::string input1;
	std::string input2;

	if (argc >= 2)
	{
		input1 = std::string(argv[1]);
		input2 = std::string(argv[2]);
	}
	else
	{
		input1 = "C:/Users/Stephan/Desktop/deformedMeshFEM.abc";
		input2 = "C:/Users/Stephan/Desktop/deformedMeshFEM.abc";
	}

	std::cout << "FILE 1: " << input1 << std::endl;
	std::cout << "FILE 2: " << input2 << std::endl;

	AbcReader reader1;
	if (!reader1.openArchive(input1))
	{
		std::cout << "Error Opening first file, aborting..." << std::endl;
		return 1;
	}

	reader1.sampleSpecific(0);

	std::vector<Eigen::Vector3f>& positions = reader1.getPositions();
	for (int i = 0; i < positions.size(); ++i)
	{
		std::cout << positions[i] << std::endl;
	}

	//AbcReader reader2;
	//if (!reader2.openArchive(input2))
	//{
	//	std::cout << "Error Opening first file, aborting..." << std::endl;
	//	return 1;
	//}




	return 0;
}
