#pragma once
#include <vector>
#include <string>
#include <fstream>
#include <iostream>

#include <Eigen\Dense>

struct PBDSolverTracker
{
	PBDSolverTracker()
	{
		filenameS = "SecondPiolaKirchoffTensor.m";
		filenamePF = "FirstPiolaKirchoffTensor.m";
		filenameF = "DeformationGradient.m";
	}

	void generateFileNames(int idx, int version)
	{
		filenameS = generateSingleFileName("SecondPiolaKirchoffTensor", idx, version);
		filenamePF = generateSingleFileName("FirstPiolaKirchoffTensor", idx, version);
		filenameF = generateSingleFileName("DeformationGradient", idx, version);
		filenameAverageDeltaXLength = generateSingleFileName("DeltaXAverageLength", idx, version);
	}

	std::string generateSingleFileName(std::string baseName, int idx, int version)
	{
		std::stringstream ss;
		ss << baseName << "_" << idx << "_" << version << ".m";
		std::string result = ss.str();
		ss.clear();
		
		return result;
	}

	std::string filenameS;
	std::string filenamePF;
	std::string filenameF;

	std::string filenameAverageDeltaXLength;
	std::vector<float> averageDeltaXLength;

	std::vector<Eigen::Matrix3f> S;
	std::vector<Eigen::Matrix3f> PF;
	std::vector<Eigen::Matrix3f> F;

	void writeMatlabFile_deltaX()
	{
		std::ofstream file;
		file.open(filenameAverageDeltaXLength);
		if (!file.is_open())
		{
			std::cout << "ERROR: Could not write [ " << filenameAverageDeltaXLength << " ]." << std::endl;
			return;
		}

		file << "averageDeltaXLength = [ ";
		for (int i = 0; i < averageDeltaXLength.size(); ++i)
		{
			file << averageDeltaXLength[i];
			if (i < averageDeltaXLength.size() - 1)
			{
				file << ", ";
			}
			else
			{
				file << " ];";
			}
		}

		file.close();
		file.clear();
	}

	void writeMatlabFile_S()
	{
		array3d_2_matlab(S, filenameS, "S");
	}

	void writeMatlabFile_F()
	{
		array3d_2_matlab(F, filenameF, "F");
	}

	void writeMatlabFile_PF()
	{
		array3d_2_matlab(PF, filenamePF, "PF");
	}

	void writeAll()
	{
		std::cout << "WRITING TRACKED SOLVER DATA..." << std::endl;
		writeMatlabFile_S();
		writeMatlabFile_PF();
		writeMatlabFile_F();
		writeMatlabFile_deltaX();

		std::cout << "-------------------" << std::endl;
	}

private:
	bool array3d_2_matlab(const std::vector<Eigen::Matrix3f>& cArray, const std::string& fileName, const std::string& arrayName)
	{
		std::ofstream file;
		file.open(fileName);
		if (!file.is_open())
		{
			std::cout << "ERROR: Could not write [ " << fileName << " ]." << std::endl;
			return false;
		}

		file << arrayName << " = zeros(3, 3, " << cArray.size() << ");" << std::endl;
		for (int d = 0; d < cArray.size(); ++d)
		{
			file << arrayName << "(:, :, " << d + 1 << ") = ";
			array2d_to_matlab(cArray[d], file);
		}

		file.close();
		file.clear();

		return true;
	}


	void array2d_to_matlab(const Eigen::Matrix3f& cArray, std::ofstream& file)
	{
		file << "[ ";
		for (int row = 0; row < 3; ++row)
		{
			for (int col = 0; col < 3; ++col)
			{
				file << cArray(row, col);

				if (col == 2)
				{
					file << ";" << std::endl;
				}
				else
				{
					file << ", ";
				}
			}

		}
		file << "];" << std::endl;
	}
};