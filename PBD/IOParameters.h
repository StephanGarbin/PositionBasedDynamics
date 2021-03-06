#pragma once

#include <string>
#include <vector>

struct IOParameters
{
	std::string nodeFile;
	std::string elementFile;
	std::string constraintFile;

	std::vector<std::string> trackerFiles;

	std::vector<std::string> movingHardConstraintsIndexFiles;
	std::vector<std::string> movingHardConstraintsLocatorNames;
	std::string movingHardConstraintsAbcArchive;

	std::string customTetAttributeFile;

	void initialiseToDefaults()
	{
		nodeFile = ("barout.node");
		elementFile = ("barout.ele");
		constraintFile = ("barLowVertexConstraints.txt");

		trackerFiles.push_back("def1_constrainedEdgeBottom.txt");
		trackerFiles.push_back("def1_constrainedEdgeTop.txt");
		trackerFiles.push_back("def1_freeEdgeTop.txt");
		trackerFiles.push_back("def1_movingPin.txt");
		//trackerFiles.push_back("def1_freeEdgeBottom.txt");

		//nodeFile = ("liver_processed.1.node");
		//elementFile = ("liver_processed.1.ele");
		//constraintFile = ("pigLiver2VertexConstraints.txt");

		//nodeFile = ("liver_processedD.1.node");
		//elementFile = ("liver_processedD.1.ele");
		//constraintFile = ("pigLiver4VertexConstraints.txt");
	}
};