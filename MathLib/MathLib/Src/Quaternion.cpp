#include "Quaternion.h"

namespace Math
{
	QXquaternion::QXquaternion() : w(1), v(0, 0, 0)
	{
	}

	QXquaternion::QXquaternion(const QXquaternion& q) : w(q.w), v(q.v.x, q.v.y, q.v.z)
	{
	}

	QXquaternion::QXquaternion(const QXquaternion&& q)
	{
		std::move(q);
	}

	QXquaternion::QXquaternion(float vw, Vec3 vQ) : w(vw), v(Vec3(vQ))
	{
	}

	QXquaternion::QXquaternion(float vw, float vx, float vy, float vz) : w(vw), v(vx, vy, vz)
	{
	}

	void QXquaternion::nullQuaternion()
	{
		w = 0;
		v.x = v.y = v.z = 0;
	}

	float QXquaternion::QuaternionLength()
	{
		return sqrt(powf(w, 2) + (v.Dot(v)));
	}

	float QXquaternion::sqrtRootQuaternion()
	{
		return dotProductQuaternion(conjugateQuaternion());
	}

	std::string QXquaternion::ToString() const
	{
		std::string quat = std::to_string(w) + ", " + v.ToString();

		return quat;
	}

	QXquaternion& QXquaternion::normalizeQuaternion()
	{
		float s = 1 / QuaternionLength();

		return multQuaternion(s);
	}

	QXquaternion& QXquaternion::conjugateQuaternion()
	{
		QXquaternion res(w, -v.x, -v.y, -v.z);
		return res;
	}

	QXquaternion& QXquaternion::inverseQuaternion()
	{
		return conjugateQuaternion().multQuaternion(1 / sqrtRootQuaternion());
	}

	void QXquaternion::negateQuaternion()
	{
		w = -w;
		v.x = -v.x;
		v.y = -v.y;
		v.z = -v.z;
	}

	QXquaternion& QXquaternion::returnNegateQuaternion()
	{
		QXquaternion res = QXquaternion();

		res.w = -w;
		res.v.x = -v.x;
		res.v.y = -v.y;
		res.v.z = -v.z;

		return res;
	}

	QXquaternion& QXquaternion::addQuaternion(QXquaternion& q)
	{
		QXquaternion res = QXquaternion();

		res.w = w + q.w;
		res.v = v + q.v;
		return res;
	}

	QXquaternion& QXquaternion::multQuaternion(float s)
	{
		QXquaternion res = QXquaternion();

		res.w = s * w;
		res.v = v * s;
		return res;
	}

	QXquaternion& QXquaternion::multQuaternion(QXquaternion& q)
	{
		QXquaternion res = QXquaternion();

		res.w = w * q.w - v.x * q.v.x - v.y * q.v.y - v.z * q.v.z;
		res.v.x = w * q.v.x + v.x * q.w + v.y * q.v.z - v.z * q.v.y;
		res.v.y = w * q.v.y - v.x * q.v.z + v.y * q.w + v.z * q.v.x;
		res.v.z = w * q.v.z + v.x * q.v.y - v.y * q.v.x + v.z * q.w;

		return res;
	}

	float QXquaternion::dotProductQuaternion(QXquaternion& q)
	{
		return (w * q.w + v.Dot(q.v));
	}

	QXquaternion& QXquaternion::slerpQuaternion(QXquaternion& q, float t)
	{
		float theta{ acos(dotProductQuaternion(q)) };

		if (theta < 0.f)
			* this = *this * -1.f;

		return ((*this * sin((1 - t) * theta) + q * sin(t * theta)) * (1 / sin(theta))).normalizeQuaternion();
	}

	QXquaternion QXquaternion::slerpQuaternion(QXquaternion q1, QXquaternion q2, float t)
	{
		float theta{ q1.dotProductQuaternion(q2) };

		if (theta < 0.f)
			q1 = q1 * -1.f;

		return ((q1 * sin((1 - t) * theta) + q2 * sin(t * theta)) * (1 / sin(theta))).normalizeQuaternion();
	}

	QXquaternion& QXquaternion::subQuaternion(QXquaternion& q)
	{
		QXquaternion res = QXquaternion();

		res.w = w - q.w;
		res.v = v - q.v;
		return res;
	}

	Mat4 QXquaternion::convertQuaternionToMat()
	{
		Mat4 res;

		res.array[0] = 1 - (2 * powf(v.y, 2)) - (2 * powf(v.z, 2));
		res.array[1] = (2 * v.x * v.y) - (2 * w * v.z);
		res.array[2] = (2 * v.x * v.z) + (2 * w * v.y);

		res.array[4] = (2 * v.x * v.y) + (2 * w * v.z);
		res.array[5] = 1 - (2 * powf(v.x, 2)) - (2 * powf(v.z, 2));
		res.array[6] = (2 * v.y * v.z) - (2 * w * v.x);

		res.array[8] = (2 * v.x * v.z) - (2 * w * v.y);
		res.array[9] = (2 * v.y * v.z) + (2 * w * v.x);
		res.array[10] = 1 - (2 * powf(v.x, 2)) - (2 * powf(v.y, 2));
		res.array[15] = 1;

		return res;
	}

	Mat4 QXquaternion::convertQuaternionToMat(QXquaternion& q)
	{
		Mat4 res;

		res.array[0] = 1 - (2 * powf(q.v.y, 2)) - (2 * powf(q.v.z, 2));
		res.array[1] = (2 * q.v.x * q.v.y) - (2 * q.w * q.v.z);
		res.array[2] = (2 * q.v.x * q.v.z) + (2 * q.w * q.v.y);

		res.array[4] = (2 * q.v.x * q.v.y) + (2 * q.w * q.v.z);
		res.array[5] = 1 - (2 * powf(q.v.x, 2)) - (2 * powf(q.v.z, 2));
		res.array[6] = (2 * q.v.y * q.v.z) - (2 * q.w * q.v.x);

		res.array[8] = (2 * q.v.x * q.v.z) - (2 * q.w * q.v.y);
		res.array[9] = (2 * q.v.y * q.v.z) + (2 * q.w * q.v.x);
		res.array[10] = 1 - (2 * powf(q.v.x, 2)) - (2 * powf(q.v.y, 2));
		res.array[15] = 1;

		return res;
	}

	QXquaternion QXquaternion::convertMatToQuaternion(Mat4 m)
	{
		float qw = sqrt(1 + m.array[0] + m.array[5] + m.array[10]) / 2.f;
		float qx = (m.array[9] - m.array[6]) / (4 * qw);
		float qy = (m.array[2] - m.array[8]) / (4 * qw);
		float qz = (m.array[4] - m.array[1]) / (4 * qw);

		return QXquaternion(qw, Vec3(qx, qy, qz));
	}

	QXquaternion QXquaternion::convertEulerAngleToQuaternion(Vec3& euler)
	{
		float c1 = cos(euler.x / 2);
		float c2 = cos(euler.y / 2);
		float c3 = cos(euler.z / 2);
		float s1 = sin(euler.x / 2);
		float s2 = sin(euler.y / 2);
		float s3 = sin(euler.z / 2);

		float qw = ((c1 * c2 * c3) - (s1 * s2 * s3));
		float qx = ((s1 * s2 * c3) + (c1 * c2 * s3));
		float qy = ((s1 * c2 * c3) + (c1 * s2 * s3));
		float qz = ((c1 * s2 * c3) - (s1 * c2 * s3));

		return QXquaternion(qw, Vec3(qx, qy, qz));
	}

	QXquaternion& QXquaternion::operator=(const QXquaternion& q)
	{
		w = q.w;
		v = q.v;

		return *this;
	}

	QXquaternion& QXquaternion::operator*(float s)
	{
		QXquaternion res = QXquaternion();

		res.w = s * w;
		res.v = v * s;
		return res;
	}

	QXquaternion& QXquaternion::operator*(const QXquaternion& q)
	{
		QXquaternion res = QXquaternion();

		res.w = w * q.w - v.x * q.v.x - v.y * q.v.y - v.z * q.v.z;
		res.v.x = w * q.v.x + v.x * q.w + v.y * q.v.z - v.z * q.v.y;
		res.v.y = w * q.v.y - v.x * q.v.z + v.y * q.w + v.z * q.v.x;
		res.v.z = w * q.v.z + v.x * q.v.y - v.y * q.v.x + v.z * q.w;

		return res;
	}

	QXquaternion& QXquaternion::operator+(const QXquaternion& q)
	{
		QXquaternion res = QXquaternion();

		res.w = w + q.w;
		res.v = v + q.v;
		return res;
	}

	QXquaternion& QXquaternion::operator-(const QXquaternion& q)
	{
		QXquaternion res = QXquaternion();

		res.w = w - q.w;
		res.v = v - q.v;
		return res;
	}
}