/*
Please note that this has been submitted before for other courseworks.
*/
#pragma once
#include <string>

#include <OpenMesh\Core\IO\MeshIO.hh>
#include "OMDefinitions.h"

void loadMesh(const std::string& fileName, meshType& mesh);

void saveMesh(const std::string& fileName, meshType& mesh);