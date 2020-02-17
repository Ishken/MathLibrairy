#pragma once

#include "Vec3.h"
#include "Mat.h"
#include "Mat4.h"

namespace Math
{
	struct	QXquaternion
	{
		#pragma region Attributes
		float	w{ 0.f };
		Vec3	v = Vec3(0, 0, 0);
		#pragma endregion
	
		#pragma region Constructors/Destructor
		QXquaternion() noexcept;
		QXquaternion(const QXquaternion& q) noexcept;
		QXquaternion(const QXquaternion&& q) noexcept;
		QXquaternion(float vw, Vec3 vQ) noexcept;
		QXquaternion(float vw, float vx, float vy, float vz) noexcept;
		~QXquaternion() = default;
		#pragma endregion

		#pragma region Functions
		QXquaternion&		addQuaternion(QXquaternion& q);
		QXquaternion&		conjugateQuaternion();
		Mat4				convertQuaternionToMat();
		float				dotProductQuaternion(QXquaternion& q);
		QXquaternion&		inverseQuaternion();
		QXquaternion&		multQuaternion(float s);
		QXquaternion&		multQuaternion(QXquaternion& q);
		void				negateQuaternion();
		QXquaternion&		normalizeQuaternion();
		void				nullQuaternion();
		float				QuaternionLength();
		QXquaternion&		returnNegateQuaternion();
		float				sqrtRootQuaternion();
		QXquaternion&		slerpQuaternion(QXquaternion& q, float t);
		QXquaternion&		subQuaternion(QXquaternion& q);
		std::string			ToString() const;
		#pragma region Static Functions
		static Mat4			convertQuaternionToMat(QXquaternion& q);
		static QXquaternion	convertMatToQuaternion(Mat4 m);
		static QXquaternion	convertEulerAngleToQuaternion(Vec3& euler);
		static QXquaternion	slerpQuaternion(QXquaternion q1, QXquaternion q2, float t);
		#pragma endregion Static Functions
		#pragma region Operator Functions
		QXquaternion&		operator=(const QXquaternion& q);
		QXquaternion&		operator*(float s);
		QXquaternion&		operator*(const QXquaternion& q);
		QXquaternion&		operator+(const QXquaternion& q);
		QXquaternion&		operator-(const QXquaternion& q);
		#pragma endregion Operator Functions
		#pragma endregion Functions


	};
}