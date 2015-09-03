#pragma once
#include <vector>
#include <memory>
#include <string>

#include <Eigen\Dense>

#include "PBDParticle.h"
#include "AbcReaderTransform.h"

class MovingHardConstraints
{
public:
	MovingHardConstraints();
	~MovingHardConstraints();

	void readConstraintFiles(const std::vector<std::string>& files);

	void readFromAbc(const std::string& fileName, const std::vector<std::string>& topBottomTransformNames);

	void initialisePositionMasses(std::vector<PBDParticle>& positions);

	void updatePositions(std::vector<PBDParticle>& positions, int currentFrame, float timeStep, int locatorIdx);

	float& getSpeed() { return m_speed; }
private:
	float m_speed;
	std::shared_ptr<AbcReaderTransform> m_reader;
	std::vector<std::vector<int>> m_constraintIndices;

	std::vector<Eigen::Vector3f> m_previousPosition;
};

