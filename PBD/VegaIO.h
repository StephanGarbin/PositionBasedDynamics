#pragma once
#include <memory>
#include <vector>

#include "PBDParticle.h"
#include "PBDTetrahedra3d.h"

class VegaIO
{
public:
	VegaIO();
	~VegaIO();

	static bool writeVegFile(const std::vector<PBDTetrahedra3d>& tetrahedra,
		const std::shared_ptr<std::vector<PBDParticle>> particles, const std::string& fileName);
};

