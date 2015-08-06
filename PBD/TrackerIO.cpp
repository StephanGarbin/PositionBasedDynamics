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
		//std::cout << std::stof(xEntries[i]) << ", " << std::stof(yEntries[i]) << std::endl;
	}

	std::cout << "Successfully read " << trackingData.size() << " frames of constraint animation from [ " << fileName << " ]. " << std::endl;

	return true;
}


Eigen::Vector3f
TrackerIO::getInterpolatedConstraintPosition(const std::vector<Eigen::Vector2f>& trackingData,
	float timeStepTracker, float timeStepSolver, float currentSystemTime)
{
	Eigen::Vector3f result;

	if (currentSystemTime == 0.0f)
	{
		return Eigen::Vector3f(trackingData[0].x(), 0.0f, trackingData[0].y());
	}

	//1. Determine previous and current frame
	float timeToReturn = currentSystemTime * (timeStepTracker / timeStepSolver);

	//2. Check this is in range
	if (timeToReturn >= timeStepTracker * (float)trackingData.size())
	{
		return Eigen::Vector3f(trackingData[trackingData.size() - 1].x(), 0.0f, trackingData[trackingData.size() - 1].y());
	}

	//3. Interpolate if necessary
	int previousFrame = std::floor(timeToReturn / timeStepTracker);
	int nextFrame = std::ceil(timeToReturn / timeStepTracker);

	float factor = (timeToReturn - ((float)previousFrame * timeStepTracker)) / timeStepTracker;

	Eigen::Vector3f prev(trackingData[previousFrame].x(), 0.0f, trackingData[previousFrame].y());
	Eigen::Vector3f next(trackingData[nextFrame].x(), 0.0f, trackingData[previousFrame].y());

	//result = factor * prev + (1.0f - factor) * next;
	result = (1.0f - factor) * prev + factor * next;

	return result;
}