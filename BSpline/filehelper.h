#pragma once
#include <vector>
#include  <vec.h>
#include <fstream>
#include <sstream>
#include <string>
#include <Eigen/Dense>
using namespace std;
using namespace Eigen;

void readdata(vector<Vector3d> &points , const char* path) {
	std::ifstream file(path);
	string line;
	points.clear();
	while (getline(file, line)) {
		istringstream iss(line);
		double x, y, z, w;
		iss >> x >> y >> z >> w;
		points.push_back(Vector3d(x, y, z));
	}
}

void readSurdata(vector<vector<Vector3d>>& points, const char* path, int &u , int &v) {
	std::ifstream file(path);
	string line;
	points.clear();
	getline(file, line);
	istringstream iss(line);
	iss >> u >> v;
	for (int i = 0; i < u; ++i) {
		for (int k = 0; k < v; ++k) {
			if(!getline(file,line))
			{
				fprintf(stderr, "data not enough\n");
				exit(1);
			}
			istringstream iss(line);
			double x, y, z, w;
			iss >> x >> y >> z >> w;
			points[i].push_back(Vector3d(x, y, z));
		}
	}
	u--;
	v--;
}