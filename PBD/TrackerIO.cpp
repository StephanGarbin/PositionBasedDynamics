#include "TrackerIO.h"

#include <fstream>
#include <iostream>
#include <algorithm>

#include <boost\algorithm\string.hpp>
#include <boost\algorithm\string\trim_all.hpp>


TrackerIO::TrackerIO()
{
}


TrackerIO::~TrackerIO()
{
}


bool
TrackerIO::readTrackerAnimationNuke(const std::string& fileName, std::vector<Eigen::Vector2f>& trackingData)
{
	std::ifstream file;
	file.open(fileName);

	if (!file.is_open())
	{
		std::cout << "ERROR: Cannot read tracking data [ " << fileName << " ]" << std::endl;
		return false;
	}

	//Read the text file
	std::string fileText;
	fileText.assign(std::istreambuf_iterator<char>(file), (std::istreambuf_iterator<char>()));

	if (fileText.empty())
	{
		std::cout << "ERROR: Tracking data file [ " << fileName << " ] is empty." << std::endl;
		return false;
	}

	//1. Split x and y
	int startOfSecond = fileText.find("{", 2);

	std::string xCoordString = fileText.substr(0, startOfSecond);
	std::string yCoordString = fileText.substr(startOfSecond, fileText.size());

	boost::algorithm::trim_all(xCoordString);
	boost::algorithm::trim_all(yCoordString);

	boost::algorithm::erase_all(xCoordString, "{");
	boost::algorithm::erase_all(xCoordString, "}");

	boost::algorithm::erase_all(yCoordString, "{");
	boost::algorithm::erase_all(yCoordString, "}");

	//2. Split based on space
	std::vector<std::string> xEntries;
	boost::algorithm::split(xEntries, xCoordString, boost::is_any_of(" "));

	std::vector<std::string> yEntries;
	boost::algorithm::split(yEntries, yCoordString, boost::is_any_of(" "));

	//3. Save data
	for (int i = 2; i < xEntries.size(); ++i)
	{
		trackingData.push_back(Eigen::Vector2f(std::stof(xEntries[i]), std::stof(yEntries[i])));
	}

	std::cout << "Successfully read " << trackingData.size() << " frames of constraint animation from [ " << fileName << " ]. " << std::endl;

	return true;
}
