#pragma once
#include <vector>
#include  <vec.h>
#include <fstream>
#include <sstream>
#include <string>
#include <Eigen/Dense>
using namespace std;
using namespace Eigen;
extern vector<Vector3d> points;

void readdata(const char* path) {
	std::ifstream file(path);
	string line;
	while (getline(file, line)) {
		istringstream iss(line);
		double x, y, z, w;
		iss >> x >> y >> z >> w;
		points.push_back(Vector3d(x, y, z));
	}
}