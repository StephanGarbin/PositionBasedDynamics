#pragma once
#include <vector>
#include <string>

#include <Eigen\Dense>

class TrackerIO
{
public:
	static bool readTrackerAnimationNuke(const std::string& fileName, std::vector<Eigen::Vector2f>& trackingData);

private:
	TrackerIO();
	~TrackerIO();
};

