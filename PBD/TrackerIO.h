#pragma once
#include <vector>
#include <string>

#include <Eigen\Dense>

class TrackerIO
{
public:
	static bool readTrackerAnimationNuke(const std::string& fileName, std::vector<Eigen::Vector2f>& trackingData);

	static Eigen::Vector3f getInterpolatedConstraintPosition(const std::vector<Eigen::Vector2f>& trackingData,
		float timeStepTracker, float timeStepSolver, float currentSystemTime);

private:
	TrackerIO();
	~TrackerIO();
};

