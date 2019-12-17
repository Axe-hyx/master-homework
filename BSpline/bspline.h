#pragma once
#include <vector>
#ifndef _MSV_VER
#include <Eigen/Dense>
#else
#include <eigen3/Eigen/Dense>
#endif
#include "filehelper.h"
using namespace std;
using namespace Eigen;

struct Curve {
	vector<float> param_t;
	vector <float> knot;
	int n, degree, order;
	vector<Vector3d> points;
	MatrixXd P, D;
  int h; // control points in apporiximation
}curve;
bool debug = 1;

int command = 2;
struct Surface{
	vector<float> param_s, param_t;
	vector <float> knot_U, knot_V;
	vector<vector<Vector3d>> points;
	vector<vector<Vector3d>> Q, P;
	MatrixXd D;
	int n, m;
	int p, q; // degree
	int order_p, order_q; // order
}surface;

enum en_t {
	UNIFORM_SPACED,
	CHORD_LENGTH
};
enum en_knot {
	k_SPACED,
	AVERAGE,
	UNIVER
};


const double Nikt(int i, int k, double t)
{
	using namespace std;
	auto& x = curve.knot;
	if (k == 1)
	{
		if (t >= x[i] && t < x[i + 1]) return 1;
		else return 0;
	}
	double d1 = 0, d2 = 0;
	if (i + k - 1 < x.size() && i < x.size())
		d1 = x[i + k - 1] - x[i];
	if (i + k < x.size() && i + 1 < x.size())
		d2 = x[i + k] - x[i + 1];
	double a, b;
	if (abs(d1) < 1e-6) a = 0;
	else a = (t - x[i]) * Nikt(i, k - 1, t) / d1;
	if (abs(d2) < 1e-6 || i + k >= x.size()) b = 0;
	else b = (x[i + k] - t) * Nikt(i + 1, k - 1, t) / d2;
	return a + b;
}

float cal_N(vector<float>& knot, int i, int k, float t) {
	if (k == 1) {
		if (t >= knot[i] && t < knot[i + 1])
		{
			return 1.0;
		}
		else
		{
			return 0.0;
		}
	}
    assert(i+k < knot.size());
	float lfrac = (knot[i + k - 1] - knot[i]);
	float rfrac = (knot[i + k] - knot[i + 1]);
	float a = 0.0f, b = 0.0f;
	if (abs(lfrac) < 1e-6 && abs(rfrac) < 1e-6) {
		return 0.0;
	}
	else if (abs(lfrac) < 1e-6) {
		a = (knot[i + k] - t) / rfrac * cal_N(knot, i + 1, k - 1, t);
	}
	else if (abs(rfrac) < 1e-6) {
		b = (t - knot[i]) / lfrac * cal_N(knot, i, k - 1, t);
	}
	else {
		a = (knot[i + k] - t) / rfrac * cal_N(knot, i + 1, k - 1, t);
		b = (t - knot[i]) / lfrac * cal_N(knot, i, k - 1, t);
	}
	return a + b;
}

const double distance(const Vector3d& v1, const Vector3d& v2) {
	return sqrt((v1(0) - v2(0)) * (v1(0) - v2(0)) + (v1(1) - v2(1)) * (v1(1) - v2(1)) + (v1(2) - v2(2)) * (v1(2) - v2(2)));
}
/*
get parm t
*/
void generatePram(const vector<Vector3d>* p, vector<float>& param_t, int n, en_t type, int a, int b) {
	param_t.clear();
	double length = 0.0;
	double ct = 0;
	vector<Vector3d> points;
	if (p != nullptr)
		points = *p;
	switch (type)
	{
	case UNIFORM_SPACED:
		for (int i = 0; i <= n;++i) {
			param_t.push_back((float)a + (float)i / (n) * (b - a));
		}
		break;
	case CHORD_LENGTH:
		for (int i = 1; i <= n;++i) {
			length += distance(points[i], points[i - 1]);
		}
		for (int i = 0; i <= n;++i) {
			if (i == 0) param_t.push_back(0);
			else {
				ct += distance(points[i], points[i - 1]);
				param_t.push_back((float)ct / length);
			}

		}
		break;
	default:
		break;
	}
	if (debug) {
		cout << "T :";
		for (size_t i = 0; i < param_t.size(); ++i) {
			cout << param_t[i] << " ";
		}
		cout << endl;
	}
}

/*
	generate knot
*/
void generateKnot(vector<float>& knot, en_knot type, int n, int degree,
	const vector<float>& param_t) {
	int order = degree + 1;
	knot.clear();
	if (type == UNIVER) {
		for (int i = 1; i <= n + 1 + order;++i) {
			knot.push_back((float)i / (n + order + 1));
		}
	}
	else if (type == k_SPACED) {
		for (int i = 0; i < order;++i) {
			knot.push_back(0);
		}
		for (int i = 1;i <= n - degree;++i) {
			knot.push_back((float)i / (n - degree + 1));
		}
		for (int i = 0; i < order;++i) {
			knot.push_back(1);
		}
	}
	else if (type == AVERAGE) {
		for (int i = 0; i < order;++i) {
			knot.push_back(0);
		}
		for (int j = 1; j <= n - degree;++j) {
			double avg = 0.0;
			for (int i = j; i <= j + degree - 1;++i) {
				avg += param_t[i];
			}
			knot.push_back(avg / degree);
		}
		for (int i = 0; i < order;++i) {
			knot.push_back(1);
		}
	}
	if (debug) {
		cout << "KNOT: (" << knot.size()<<"):";
		for (size_t i = 0; i < knot.size(); ++i) {
			cout << knot[i] << " ";
		}
		cout << endl;
	}
}

void drawPoly(int i, int k, float t) {
	float step = 0.1f;
	fprintf(stderr, "N(%d,%d): ", i, k);
	for (float c = 0; c < 1;c += step) {
		fprintf(stderr, "%f ", cal_N(curve.knot, i, k, c));
	}
	fprintf(stderr, "%f ", cal_N(curve.knot, i, k, t));
	fprintf(stderr, "\n");
}
void calColpoints(vector<Vector3d>& points, int order, vector<float >& knot, vector<float >& param, int n, MatrixXd& P) {
	MatrixXd N(n + 1, n + 1);
	for (int i = 0; i <= n; ++i) {
		for (int j = 0; j <= n; ++j) {
			N(i, j) = cal_N(knot, j, order, param[i]);
		}
	}
	N(n, n) = 1;
	//cout << "N: " << endl;
	//cout << N << endl;
	MatrixXd D(n + 1, 3);
	for (int i = 0; i <= n; ++i) {
		for (int j = 0; j < 3; ++j) {
			D(i, j) = points[i](j);
		}
	}
	//cout << "D: " << endl;
	//cout << D << endl;
	P.resize(n + 1, 3);
	P = N.inverse() * D;
	//  cout << "P:1 " << endl;
	//cout << P << endl;
	curve.D = D;
	curve.P = P;
}
void evaluate(Curve& cv) {
	MatrixXd N(cv.n + 1, cv.n + 1);
	for (int i = 0; i <= cv.n; ++i) {
		for (int j = 0; j <= cv.n; ++j) {
			N(i, j) = cal_N(cv.knot, j, cv.order, cv.param_t[i]);
		}
	}
	N(cv.n, cv.n) = 1;
	cout << "N: " << endl;
	cout << N << endl;
	cv.D.resize(cv.n + 1, 3);
	for (int i = 0; i <= cv.n; ++i) {
		for (int j = 0; j < 3; ++j) {
			cv.D(i, j) = cv.points[i](j);
		}
	}
	cout << "D: " << endl;
	cout << cv.D << endl;
	cv.P.resize(cv.n + 1, 3);
	cv.P = N.inverse() * cv.D;
	cout << "P: " << endl;
	cout << cv.P << endl;
}
void surface_approximation(Curve &cv) {
  
}
void curve_approximation(Curve &cv) {
  MatrixXd Pd(cv.h-1, 3);
  MatrixXd Qk(cv.n-1, 3);
  MatrixXd D(cv.n+1, 3);
  for(int i = 0; i < cv.n + 1; ++i) {
    for(int k = 0; k < 3; ++k) {
      D(i,k) = cv.points[i](k);
    }
  }
  cout<<"D: "<<endl;
  cout<<D<<endl;
  MatrixXd P(cv.h + 1, 3);
  for (int k = 0; k < 3; ++k) {
    P(0, k) = D(0,k);
    P(cv.h, k) = D(cv.n, k);
  }
  cout<<"P："<<endl;
  cout<<P<<endl;
  cout<<cv.n<<" "<<cv.h<<endl;
  for (int k = 1; k < cv.n; ++k) {
    for(int i = 0; i<3; ++i){
      Qk(k-1,i) = D(k,i)
      - cal_N(cv.knot, 0, cv.order, cv.param_t[k]) * P(0,i)
          - cal_N(cv.knot, cv.h, cv.order, cv.param_t[k]) * P(cv.h, i);
    }
  }
  cout << "Qk：" << endl;
  cout<<Qk<<endl;
  MatrixXd Q(cv.h - 1, 3);
  for (int i = 0; i < cv.h-1; ++i) {
    for (int k = 0; k < 3; ++k) {
      double sum = 0.0f;
      for (int j = 1; j < cv.n; ++j) {
        sum += cal_N(cv.knot, i, cv.order, cv.param_t[j]) * Qk(j-1, k);
      }
      Q(i,k) = sum;
    }
  }
  cout<<"Q: "<<endl;
  cout<<Q<<endl;
  MatrixXd N(cv.n - 1, cv.h - 1);
  for (int k = 1; k < cv.n; ++k) {
    for (int i = 0; i < cv.h-1; ++i) {
      N(k-1,i) = cal_N(cv.knot, i+1, cv.order, cv.param_t[k]);
    }
  }
  cout << "N: " << endl;
  cout << N << endl;
  MatrixXd M = N.transpose()*N;
  cout << "M: " << endl;
  cout << M<< endl;
  Pd = M.inverse() * Q;
  for (int i = 1; i < cv.h; ++i) {
    for (int k = 0; k < 3; ++k) {
      P(i, k ) = Pd(i-1, k);
    }
  }
  cout << "P: " << endl;
  cout << P << endl;
  cv.D = D;
  cv.P = P;
}
void global_curve_approximation() {
  command = 2;
  readdata(curve.points,
           "/home/switch/work/master-homework/BSpline/Release/curve_data.in");
  curve.n = curve.points.size() - 1;
  curve.degree = 2; // degree p
  curve.order = curve.degree + 1;
  curve.h = curve.n-1;
  if (curve.h >= curve.n) {
    fprintf(stderr, "ERROR: n > h\n");
    exit(1);
  }
  if (curve.h < curve.degree) {
    fprintf(stderr, "h%d %dERROR: h >= p\n", curve.h, curve.degree);
    exit(2);
  }
  generatePram(&curve.points, curve.param_t, curve.n, en_t::UNIFORM_SPACED, 0, 1);
  generateKnot(curve.knot, en_knot::k_SPACED, curve.h, curve.degree, curve.param_t); // knot size h + 1
  curve_approximation(curve);
}
void global_curve_interpolation() {
  command = 1;
  readdata(curve.points,
           "/home/switch/work/master-homework/BSpline/Release/curve_data.in");
  curve.n = curve.points.size() - 1;
  curve.degree = 3; // degree p;
  curve.order = curve.degree + 1;
  if (curve.n < curve.degree) {
    fprintf(stderr, "need degree +1 control points\n");
    exit(1);
	}
	generatePram(&curve.points, curve.param_t, curve.n, en_t::UNIFORM_SPACED, 0, 1);
	cout << "n+1: " << curve.n + 1 << " degree " << curve.degree << endl;
	generateKnot(curve.knot, en_knot::k_SPACED, curve.n, curve.degree,
		curve.param_t);
	evaluate(curve);
}

void inter_Q() {
	auto& s = surface;
	int row = surface.m, col = surface.n;
	s.Q.resize(col + 1);
	MatrixXd D(row + 1, 3);
	for (int i = 0; i <= col; ++i) {
		vector<Vector3d> points;
		for (int k = 0; k <= row; ++k) {
			points.push_back(s.points[k][i]);
			for (int c = 0; c < 3; ++c) {
				D(k, c) = s.points[k][i](c);
			}
		}
		MatrixXd Q(row + 1, 3);
		calColpoints(points, s.order_p, s.knot_U, s.param_s, row, Q); // for m;
		//calColpoints(points, s.order_q, s.knot_V, s.param_t, row, Q); // for n;
		//cout << Q << endl;
		//cout << endl;
		for (int k = 0; k <= row; ++k) {
			s.Q[i].push_back(Vector3d(Q(k, 0), Q(k, 1), Q(k, 2)));
		}
	}

	s.P.resize(row + 1);
	for (int c = 0; c <= row; ++c) {
		vector<Vector3d> points;
		MatrixXd Q(col + 1, 3);
		for (int k = 0;k <= col; ++k) {
			auto p = s.Q[k][c];
			points.push_back(p);
			for (int i = 0; i < 3;++i) {
				Q(k, i) = p(i);
			}
		}
		MatrixXd P(col + 1, 3);
		calColpoints(points, s.order_q, s.knot_V, s.param_t, col, P);
		cout << P << endl;
		cout << endl;
		for (int k = 0; k <= col; ++k) {
			s.P[c].push_back(Vector3d(P(k, 0), P(k, 1), P(k, 2)));
		}
	}
}
void globel_surface_interpolation() {
  readSurdata(surface.points,
              "/home/switch/work/master-homework/BSpline/Release/surface_data.in",
              surface.m, surface.n);
  surface.p = 2;
  surface.q = 2;
  surface.order_p = surface.p + 1;
  surface.order_q = surface.q + 1;
  if (surface.m < surface.p || surface.n < surface.q) {
    fprintf(stderr, "need degree +1 control points\n");
    exit(1);
	}
	generatePram(nullptr, surface.param_s, surface.m, en_t::UNIFORM_SPACED, 0, 1);
	generateKnot(surface.knot_U, en_knot::k_SPACED, surface.m,
		surface.p, surface.param_s);


	generatePram(nullptr, surface.param_t, surface.n, en_t::UNIFORM_SPACED, 0, 1);
	generateKnot(surface.knot_V, en_knot::k_SPACED, surface.n,
		surface.q, surface.param_t);
	inter_Q();

}

void testcase() {
  //global_curve_interpolation();
	//globel_surface_interpolation();
  global_curve_approximation();
}
void drawsurB(Surface& s) {
	glPointSize(4);
	glColor3f(0, 0, 1);
	glBegin(GL_POINTS);
	float step = 0.01f;
	for (float u = 0.0; u <= 1.0; u += step) {
		for (float v = 0.0; v <= 1.0; v += step) {
			Vector3d ve(0, 0, 0);
			for (int i = 0; i <= s.m ; ++i) {
				for (int k = 0; k <= s.n ; ++k) {
					ve += Vector3d(s.P[i][k]) * cal_N(s.knot_U, i, s.order_p, u) 
						* cal_N(s.knot_V, k, s.order_q, v);
				}
			}
			glVertex3f(ve(0), ve(1), ve(2));
		}
	}
	glEnd();
	glPointSize(10);
	glColor3f(1, 0, 0);
	glBegin(GL_POINTS);
	for (int  i = 0; i <= s.m; i++) {
		for (int k = 0; k <= s.n; ++k) {
			glVertex3f(s.P[i][k](0), s.P[i][k](1), s.P[i][k](2));
		}
	}
	glEnd();
	glPointSize(10);
	glColor3f(0, 1, 0);
	glBegin(GL_POINTS);
	for (int  i = 0; i <= s.m; i++) {
		for (int k = 0; k <= s.n; ++k) {
			glVertex3f(s.points[i][k](0), s.points[i][k](1), s.points[i][k](2));
		}
	}
	glEnd();
}
void drawcurB_app(Curve &cv) {
    glPointSize(4);
    glColor3f(0,0,1);
    glBegin(GL_POINTS);
    float step = 0.01f;
    for (float t = 0.0; t <= 1.0f; t += step) {
      Vector3d v(0,0,0);
      for (int i = 0; i <= cv.h; ++i) {
        if(i!=2 && i!= 0) continue;
        v+= Vector3d(cv.P(i,0), cv.P(i,1), cv.P(i,2)) *
          cal_N(cv.knot, i, cv.order, t);
      }
      glVertex3f(v(0), v(1), v(2));
    }
    glEnd();

    glPointSize(10);
    glColor3f(1, 0, 0);
    glBegin(GL_POINTS);
    for (int t = 0; t <= cv.n; t++) {
      glVertex3f(cv.D(t, 0), cv.D(t, 1), cv.D(t, 2));
    }
    glEnd();

    glPointSize(10);
    glColor3f(0, 1, 0);
    glBegin(GL_POINTS);
    for (int t = 0; t <= cv.h; t++) {
      glVertex3f(cv.P(t, 0), cv.P(t, 1), cv.P(t, 2));
    }
    glEnd();
}
void drawcuvB(Curve& cv) {
	glPointSize(4);
	glColor3f(0, 0, 1);
	glBegin(GL_POINTS);
	float step = 0.01f;
	for (float t = 0.0; t <= 1.0; t += step) {
		Vector3d v(0, 0, 0);
		for (int i = 0; i < cv.n + 1; ++i) {
			v += Vector3d(cv.P(i, 0), cv.P(i, 1), cv.P(i, 2)) * cal_N(cv.knot, i, cv.order, t);
		}
		glVertex3f(v(0), v(1), v(2));
	}
	glEnd();

	glPointSize(10);
	glColor3f(1, 0, 0);
	glBegin(GL_POINTS);
	for (int  t = 0; t <= cv.n; t++) {
		glVertex3f(cv.D(t, 0), cv.D(t, 1), cv.D(t, 2));
	}
	glEnd();
	glPointSize(10);
	glColor3f(0, 1, 0);
	glBegin(GL_POINTS);
	for (int t = 0; t <= cv.n; t++) {
		glVertex3f(cv.P(t, 0), cv.P(t, 1), cv.P(t, 2));
	}
	glEnd();
}

void displayB() {
	auto cv = curve;
	glColor3f(1, 1, 1);
	glBegin(GL_LINES);
	glVertex3d(0, 0, 0);
	glVertex3d(100, 0, 0);
	glVertex3d(0, 0, 0);
	glVertex3d(0, 100, 0);
	glVertex3d(0, 0, 0);
	glVertex3d(0, 0, 100);
	glEnd();

	glColor3f(1, 0, 1);
	if (curve.points.size() > 0) {
		glBegin(GL_LINES);
		for (size_t i = 1; i < cv.points.size(); ++i) {
			glVertex3d(cv.points[i - 1](0), cv.points[i - 1](1), cv.points[i - 1](2));
			glVertex3d(cv.points[i](0), cv.points[i](1), cv.points[i](2));
		}
		glEnd();
    if(command ==1)
      drawcuvB(cv);
    else
      drawcurB_app(cv);
	}
	else {
		for (int k = 0; k <= surface.m; ++k) {
			glColor3f(0.2f, 0.2f, 0.2f);
			glBegin(GL_LINES);
			auto points = surface.points[k];
			for (int i = 1; i <= surface.n;++i) {
				glVertex3d(points[i - 1](0), points[i - 1](1), points[i - 1](2));
				glVertex3d(points[i](0), points[i](1), points[i](2));
			}
			glEnd();
		}

		for (int i = 0; i <= surface.n;++i) {
			glColor3f(0.2f, 0.2f, 0.2f);
			glBegin(GL_LINES);
			for (int k = 1; k <= surface.m; ++k) {
			auto p1 = surface.points[k];
			auto p2 = surface.points[k - 1];
			glVertex3d(p2[i](0), p2[i](1), p2[i](2));
			glVertex3d(p1[i](0), p1[i](1), p1[i](2));
			}
			glEnd();
		}

		auto m = surface.m;
		auto n = surface.n;

		glPointSize(4);
		glColor3f(1, 0, 0);
		glBegin(GL_POINTS);
		auto points = surface.Q;
		for (auto pl : points) {
			for (auto p : pl) {
				glVertex3d(p(0), p(1), p(2));
			}
		}
		glEnd();
		drawsurB(surface);
	}
	glFlush();
}
