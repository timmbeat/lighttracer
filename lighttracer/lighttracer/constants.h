/**
* @file   constants.h
* @author Sebastian Maisch <sebastian.maisch@googlemail.com>
* @date   2017.10.04
*
* @brief  Definitions for several constants used.
*/

#pragma once

#include <algorithm>

namespace slabProfiles
{

	using Real = long double;

	template<typename Real> constexpr Real pi()
	{
		return static_cast<Real>(3.14159265358979323846);
	}
	template<typename Real> constexpr Real one_over_four_pi()
	{
		return static_cast<Real>(0.07957747154594766788444188168626);
	}
	template<typename Real> constexpr Real cos90_d()
	{
		return static_cast<Real>(1 / 1000000000);
	}
	template<typename Real> constexpr Real cos_1()
	{
		return static_cast<Real>(0.999999999999);
	}
	template<typename Real> constexpr Real clamp(Real x, Real x0, Real x1)
	{
		return std::min(std::max(x0, x), x1);
	}

	template<typename Real> constexpr Real lerp(Real v0, Real v1, Real t)
	{
		return v0 * (1.0 - t) + v1 * t;
	}

	namespace corrections
	{

		template<typename Real>  Real CalcR(Real eta)
		{
			Real meta = 1.0 - eta;
			Real peta = 1.0 + eta;
			Real diveta = meta / peta;
			return diveta * diveta;
		}

		template<typename Real>  std::pair<Real, Real> CalcSourceMod(Real depth, Real R, Real rhop, Real sigmap_t)
		{
			Real origSource = rhop * sigmap_t * std::exp(-sigmap_t * depth);
			Real RS = R * origSource;
			Real Q1 = 1.0 / (1.0 - RS * RS);
			Real Q2 = RS * Q1;
			return std::pair<Real, Real>(Q1, Q2);
		}
	}
} // namespace deepProfile
