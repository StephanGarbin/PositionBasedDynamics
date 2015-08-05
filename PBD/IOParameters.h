#pragma once

#include <string>

struct IOParameters
{
	std::string nodeFile;
	std::string elementFile;
	std::string constraintFile;

	void initializeToDefaults()
	{
		nodeFile = ("barout.node");
		elementFile = ("barout.ele");
		constraintFile = ("barLowVertexConstraints.txt");

		//nodeFile = ("liver_processed.1.node");
		//elementFile = ("liver_processed.1.ele");
		//constraintFile = ("pigLiver2VertexConstraints.txt");

		//nodeFile = ("liver_processedD.1.node");
		//elementFile = ("liver_processedD.1.ele");
		//constraintFile = ("pigLiver4VertexConstraints.txt");
	}
};