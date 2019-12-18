#pragma once
#include <vector>
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
		iss >> x >> y >> z;
		points.push_back(Vector3d(x, y, z));
	}
}

void readSurdata(vector<vector<Vector3d>>& points, const char* path, int &u , int &v) {
	std::ifstream file(path);
	string head;
	points.clear();
	getline(file, head);
	istringstream iss(head);
	iss >> u >> v;
	points.resize(u);
	for (int i = 0; i < (int)u; ++i) {
		for (int k = 0; k < (int)v; ++k) {
			string line;
			if(!getline(file,line))
			{
				fprintf(stderr, "data not enough\n");
				exit(1);
			}
			istringstream is(line);
			double x, y, z;
			is >> x >> y >> z ;
			points[i].push_back(Vector3d(x, y, z));
		}
	}
	u--;
	v--;
}
