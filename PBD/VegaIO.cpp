#include "VegaIO.h"
#include <fstream>
#include <iostream>


VegaIO::VegaIO()
{
}


VegaIO::~VegaIO()
{
}

bool
VegaIO::writeVegFile(const std::vector<PBDTetrahedra3d>& tetrahedra,
const std::shared_ptr<std::vector<PBDParticle>> particles, const std::string& fileName)
{
	std::ofstream file;
	file.open(fileName);

	if (!file.is_open())
	{
		std::cout << "ERROR: Could not write file: " << fileName << std::endl;
	}

	//1. Write Particles
	file << "*VERTICES" << std::endl;
	file << particles->size() << " 3" << std::endl;
	for (int i = 0; i < particles->size(); ++i)
	{
		file << i << " " << (*particles)[i].position()[0]
			<< " " << (*particles)[i].position()[1]
			<< " " << (*particles)[i].position()[2] << std::endl;
	}

	//2. Write Tetrahedra
	file << "*ELEMENTS" << std::endl;
	file << "TETS" << std::endl;
	file << tetrahedra.size() << " 4" << std::endl;
	for (int i = 0; i < tetrahedra.size(); ++i)
	{
		file << i << " " << tetrahedra[i].getVertexIndices()[0]
			<< " " << tetrahedra[i].getVertexIndices()[1]
			<< " " << tetrahedra[i].getVertexIndices()[2]
			<< " " << tetrahedra[i].getVertexIndices()[3] << std::endl;
	}

	file.close();

	std::cout << "Successfully written: [" << fileName << "] to disk." << std::endl;
	return true;
} 
