#pragma once
#include <vector>
#include <eigen3/Eigen/Dense>
#include "filehelper.h"
using namespace std;
using namespace Eigen;

struct Curve {
	vector<float> param_t;
	vector <float> knot;
	int n, degree, order;
	vector<Vector3d> points;
	MatrixXd P, D;
}curve;
bool debug = 1;

struct {
	vector<float> param_s, param_t;
	vector <float> knot_U, knot_V;
	vector<vector<Vector3d>> points;
  vector<vector<Vector3d>> Q;
	MatrixXd P, D;
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
  if(p!=nullptr)
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
void generateKnot(vector<float>&knot, en_knot type, int n, int degree, 
	const vector<float> &param_t) {
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
		cout << "KNOT: " << knot.size();
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
void calColpoints(vector<Vector3d> &points, int order, vector<float > & knot, vector<float > & param, int n, MatrixXd &P) {
  MatrixXd N(n + 1, n + 1);
  for (int i = 0; i <= n; ++i) {
    for (int j = 0; j <= n; ++j) {
      N(i, j) = cal_N(knot, j, order, param[i]);
    }
  }
  N(n, n) = 1;
  //cout << "N: " << endl;
  //cout << N << endl;
  MatrixXd D(n+1, 3);
  for (int i = 0; i <= n; ++i) {
    for (int j = 0; j < 3; ++j) {
      D(i, j) = points[i](j);
    }
  }
  //cout << "D: " << endl;
  //cout << D << endl;
  P.resize(n+1, 3);
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
void global_curve_interpolation() {
	readdata(curve.points, "D:\\working\\master-homework\\BSpline\\Release\\curve_data.in");
	curve.n = curve.points.size() - 1;
	curve.degree = 3;
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
	auto s = surface;
  int m = s.m;
  cout<<m;
  MatrixXd N(m + 1, m + 1);
  for (int i = 0; i <= m; ++i) {
    for (int k = 0; k <= m; ++k) {
      N(i, k) = cal_N(s.knot_U, k, s.order_p, s.param_s[i]);
    }
  }
  N(s.m, s.m) = 1;
  cout << "N is " << N << endl;
  const int n = s.n;
  s.Q.resize(n + 1);
  MatrixXd D(m + 1, 3);
  for (int i = 0; i <= n; ++i) {
    vector<Vector3d> points;
    for (int k = 0; k <= m; ++k) {
      points.push_back(s.points[k][i]);
      for (int c = 0; c < 3; ++c) {
        D(k, c) = s.points[k][i](c);
      }
    }
    MatrixXd Q(m + 1, 3);
    calColpoints(points, s.order_p, s.knot_U, s.param_s, m, Q);
    cout << Q << endl;
    cout<<endl;
    for (int k = 0; k <= m; ++k) {
      s.Q[i].push_back(Vector3d(Q(k,0), Q(k,1), Q(k,2)));
    }
  }
}
void globel_surface_interpolation() {
  readSurdata(surface.points,
              "/home/switch/work/master-homework/BSpline/Release/surface_data.in",
              surface.m, surface.n);
  surface.p = 2;
  surface.q = 2;
  surface.order_p = surface.p+1;
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
	globel_surface_interpolation();
}
void drawBspline(Curve& cv) {
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
	for (float t = 0; t <= cv.n; t++) {
		glVertex3f(cv.D(t, 0), cv.D(t, 1), cv.D(t, 2));
	}
	glEnd();
	glPointSize(10);
	glColor3f(0, 1, 0);
	glBegin(GL_POINTS);
	for (float t = 0; t <= cv.n; t++) {
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
		drawBspline(cv);
	}
	else {
		for (int k = 0; k <= surface.m; ++k) {
			glBegin(GL_LINES);
			auto points = surface.points[k];
			for (int i = 1; i <= surface.n;++i) {
				glVertex3d(points[i - 1](0), points[i - 1](1), points[i - 1](2));
				glVertex3d(points[i](0), points[i](1), points[i](2));
			}
			glEnd();
		}
    auto m = surface.m;
    auto n = surface.n;

    glPointSize(4);
    glBegin(GL_POINTS);
    auto points = surface.Q;
    for (int i = 0; i <= m; ++i) {
      for (int k = 0; k <= n; ++k) {
        //glVertex3d(points[i][k](0), points[i][k](1), points[i][k](2));

      }
    }
    glEnd();
        }


	glFlush();
}
