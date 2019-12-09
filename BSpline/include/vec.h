#pragma once

#include <cmath>
#include <ostream>

class vec3f {
public:
	union {
		struct {
			float _x, _y, _z;
		};
		struct {
			float _v[3];
		};
	};

	vec3f() {
		_x = 0.0;
		_y = 0.0;
		_z = 0.0;
	}
	vec3f(float x, float y, float z) {
		_x = x;
		_y = y;
		_z = z;
	}
	vec3f(const vec3f& v) {
		_x = v._x;
		_y = v._y;
		_z = v._z;
	}
	const float x() const { return _x; }
	const float y() const { return _y; }
	const float z() const { return _z; }
	bool operator ==(const vec3f& v) {
		return (_x == v._x && _y == v._y && _z == v._z);
	}
	friend std::ostream& operator<<(std::ostream& out, const vec3f& v) {
		out << "v: (" << v._x << " " << v._y << " " << v._z << ")" << std::endl;
		return out;
	}
	vec3f operator*(float t) { return vec3f(_x * t, _y * t, _z * t); }
	vec3f operator-(const vec3f& v2) {
		return vec3f(_x - v2._x, _y - v2._y, _z - v2._z);
	}
	vec3f operator-() {
		return vec3f(-_x, -_y, -_z);
	}
	vec3f operator+(float t) { return vec3f(_x + t, _y + t, _z + t); }

	float length() { return _x * _x + _y * _y + _z * _z; }

	vec3f& unit() {
		float t = 1 / sqrt(length());
		_x *= t;
		_y *= t;
		_z *= t;
		return *this;
	}
	friend vec3f cross(const vec3f& v1, const vec3f& v2) {
		return vec3f(v1._y * v2._z - v1._z * v2._y,
			-(v1._x * v2._z - v1._z * v2._x),
			v1._x * v2._y - v1._y * v2._x);
	}
	friend float dot(const vec3f& v1, const vec3f& v2) {
		return v1._x * v2._x + v1._y * v2._y +
			v1._z * v2._z;
	}
};
