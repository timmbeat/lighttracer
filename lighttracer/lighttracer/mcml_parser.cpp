/**
* @file   mcml_parser.cpp
* @author Sebastian Maisch <sebastian.maisch@uni-ulm.de>
* @date   2018.02.07
*
* @brief  Parser for MCML output files.
*/

#include "stdafx.h"
#include "mcml_parser.h"
#include "constants.h"
#include <fstream>
#include <sstream>
#include <cassert>
#include <iostream>

namespace mcml
{

	MCMLParser::MCMLParser(const std::string & filename) :
		filename_{ filename }
	{
		std::ifstream mcoIn{ filename_ };
		std::string line;

		std::getline(mcoIn, line);
		if (line.substr(0, 2) != "A1") throw std::runtime_error("Wrong file format.");

		while (std::getline(mcoIn, line))
		{
			std::istringstream iss(line);
			std::string blockName;
			iss >> blockName;

			if (blockName == "InParm" || blockName == "InParam")
				ProcessInParams(mcoIn);
			else if (blockName == "RAT")
				ProcessRAT(mcoIn);
			else if (blockName == "A_l")
				ProcessFloatList(mcoIn, absorptionsLayer_);
			else if (blockName == "A_z")
				ProcessFloatList(mcoIn, absorptionsZ_, 0.1);
			else if (blockName == "Rd_r")
				ProcessFloatList(mcoIn, rProfileR_, 0.01);
			else if (blockName == "Rd_a")
				ProcessFloatList(mcoIn, rProfileA_);
			else if (blockName == "Tt_r")
				ProcessFloatList(mcoIn, tProfileR_, 0.01);
			else if (blockName == "Tt_a")
				ProcessFloatList(mcoIn, tProfileA_);
			else if (blockName == "A_rz")
				ProcessFloatList2(mcoIn, absorptionsRZ_, 0.001);
			else if (blockName == "Rd_ra")
				ProcessFloatList2(mcoIn, rProfileRA_, 0.01);
			else if (blockName == "Tt_ra")
				ProcessFloatList2(mcoIn, tProfileRA_, 0.01);
		}
	}

	std::size_t MCMLParser::get_numr() const
	{
		return numBinsR_;
	}

	std::size_t MCMLParser::get_numa() const
	{
		return numBinsA_;
	}

	std::size_t MCMLParser::get_numphotons() const
	{
		return numPhotons_;
	}

	double MCMLParser::get_dr_() const
	{
		return dr_;
	}

	double MCMLParser::get_dz_() const
	{
		return dz_;
	}

	double MCMLParser::get_diffuse_reflectance() const
	{
		return diffuseReflectance_;
	}

	void MCMLParser::SmoothData()
	{
		rProfileRSmooth_.resize(rProfileR_.size());
		tProfileRSmooth_.resize(tProfileR_.size());

		assert(rProfileR_.size() == tProfileR_.size());

		std::size_t width = 2;
		std::size_t part = 1;
		bool grow = true;
		for (std::size_t i = 0; i < rProfileR_.size(); ++i)
		{
			if (i <= 10)
			{
				rProfileRSmooth_[i] = rProfileR_[i];;
				tProfileRSmooth_[i] = tProfileR_[i];
				continue;
			}
			double weightSum = 0.0;
			for (std::size_t j = i - width; j <= i + width; ++j)
			{
				double weight = 1.0;
				if (j >= rProfileR_.size() - 1)
				{
					rProfileRSmooth_[i] += weight * rProfileR_[rProfileR_.size() - 2];
					tProfileRSmooth_[i] += weight * tProfileR_[rProfileR_.size() - 2];
				}
				else
				{
					rProfileRSmooth_[i] += weight * rProfileR_[j];
					tProfileRSmooth_[i] += weight * tProfileR_[j];
				}

				weightSum += weight;
			}

			rProfileRSmooth_[i] /= weightSum;
			tProfileRSmooth_[i] /= weightSum;

			if (grow && i > part * rProfileR_.size() / (part + 1))
			{
				part += 1;
				width *= 2;

				if (width >= 32) grow = false;
			}
		}
	}

	void MCMLParser::ProcessInParams(std::ifstream& mcoIn)
	{
		std::string testFilename;

		std::string line;
		std::getline(mcoIn, line);
		std::istringstream issFilename(line);
		issFilename >> testFilename;
		// assert(testFilename == filename_);

		std::getline(mcoIn, line);
		std::istringstream issPhotons(line);
		issPhotons >> numPhotons_;

		std::getline(mcoIn, line);
		std::istringstream issBinSizes(line);
		issBinSizes >> dz_ >> dr_;
		dz_ *= 10.0;
		dr_ *= 10.0;

		std::getline(mcoIn, line);
		std::istringstream issNumBins(line);
		issNumBins >> numBinsZ_ >> numBinsR_ >> numBinsA_;

		zValues_.resize(numBinsZ_);
		absorptionsZ_.resize(numBinsZ_);
		for (std::size_t i = 0; i < zValues_.size(); ++i) zValues_[i] = (static_cast<double>(i) + 0.5) * dz_;

		rValues_.resize(numBinsR_);
		rProfileR_.resize(numBinsR_);
		tProfileR_.resize(numBinsR_);
		absorptionsRZ_.resize(numBinsR_);
		rProfileRA_.resize(numBinsR_);
		tProfileRA_.resize(numBinsR_);
		for (std::size_t i = 0; i < rValues_.size(); ++i)
		{
			rValues_[i] = (static_cast<double>(i) + 0.5) * dr_;
			absorptionsRZ_[i].resize(numBinsZ_);
			rProfileRA_[i].resize(numBinsA_);
			tProfileRA_[i].resize(numBinsA_);
		}

		aValues_.resize(numBinsA_);
		rProfileA_.resize(numBinsA_);
		tProfileA_.resize(numBinsA_);
		double da = (slabProfiles::pi<double>() * 0.5) / static_cast<double>(numBinsA_);
		for (std::size_t i = 0; i < aValues_.size(); ++i) aValues_[i] = (static_cast<double>(i) + 0.5) * da;

		// skip empty parts until layers start.
		line = SkipEmptyLines(mcoIn);
		std::istringstream issNumLayers(line);
		std::size_t numLayers;
		issNumLayers >> numLayers;
		layers_.resize(numLayers);
		absorptionsLayer_.resize(numLayers);

		line = SkipEmptyLines(mcoIn);
		std::istringstream etaTLine(line);
		etaTLine >> etaTop_;

		for (auto& layer : layers_)
		{
			std::getline(mcoIn, line);
			std::istringstream layerLine(line);
			layerLine >> layer.eta_ >> layer.mua_ >> layer.mus_ >> layer.g_ >> layer.depth_;
			layer.mua_ /= 10.0;
			layer.mus_ /= 10.0;
			layer.depth_ *= 10.0;
		}

		std::getline(mcoIn, line);
		std::istringstream etaBLine(line);
		etaBLine >> etaBelow_;
	}

	void MCMLParser::ProcessRAT(std::ifstream& mcoIn)
	{
		std::string line;
		std::getline(mcoIn, line);
		std::istringstream srLine(line);
		srLine >> specularReflectance_;

		std::getline(mcoIn, line);
		std::istringstream drLine(line);
		drLine >> diffuseReflectance_;

		std::getline(mcoIn, line);
		std::istringstream afLine(line);
		afLine >> absorbedFract_;

		std::getline(mcoIn, line);
		std::istringstream ttLine(line);
		ttLine >> transmittance_;
	}

	void MCMLParser::ProcessFloatList(std::ifstream& mcoIn, std::vector<Real>& floats, Real factor)
	{
		std::string line;
		std::size_t i = 0;

		while (std::getline(mcoIn, line))
		{
			if (line.empty()) break;

			std::istringstream iss(line);
			double f;
			iss >> f;
			floats[i] = f * factor;
			i += 1;
		}

		assert(i == floats.size());
	}

	void MCMLParser::ProcessFloatList2(std::ifstream & mcoIn, std::vector<std::vector<Real>>& floats, Real factor)
	{
		std::string line;
		std::size_t i = 0, j = 0;

		while (std::getline(mcoIn, line))
		{
			if (line.empty() || line.substr(0, 1) == "#") break;

			std::istringstream iss(line);
			while (iss)
			{
				double f;
				iss >> f;
				if (iss)
				{
					floats[i][j] = f * factor;
					NextList2(i, j, floats);
				}
			}
		}

		assert(i == floats.size() && j == 0);
	}

	void MCMLParser::NextList2(std::size_t& i, std::size_t& j, const std::vector<std::vector<Real>>& floats)
	{
		j += 1;
		if (floats[i].size() <= j)
		{
			j = 0;
			i += 1;
		}
	}

	std::vector<Real> MCMLParser::get_ra() const
	{
		return rProfileR_;
	}

	std::string MCMLParser::SkipEmptyLines(std::ifstream& mcoIn)
	{
		std::string line;
		std::getline(mcoIn, line);
		while (line.empty() || line.substr(0, 1) == "#") std::getline(mcoIn, line);
		return line;
	}
}
//	template<slabProfiles::ProfileProperties PROPERTIES> MCMLProfile<PROPERTIES>::MCMLProfile(Real alphap, Real sigmapt, Real g, Real eta, Real depth,
//																							  const std::vector<Real>& radiusSamples,
//																							  const std::vector<Real>& rValues, const std::vector<Real>& tValues) :
//		DepthProfileGeneric{ static_cast<Real>(alphap), static_cast<Real>(sigmapt),
//		static_cast<Real>(g), static_cast<Real>(eta), static_cast<Real>(depth), radiusSamplesInternal_, false }
//	{
//		radiusSamplesInternal_.resize(radiusSamples.size());
//		ProfileValues<slabProfiles::ProfileType::R>().clear();
//		ProfileValues<slabProfiles::ProfileType::R>().resize(radiusSamplesInternal_.size());
//		ProfileValues<slabProfiles::ProfileType::T>().clear();
//		ProfileValues<slabProfiles::ProfileType::T>().resize(radiusSamplesInternal_.size());
//
//		for (std::size_t i = 0; i < radiusSamples.size(); ++i)
//		{
//			radiusSamplesInternal_[i] = static_cast<Real>(radiusSamples[i]);
//			Real factor = 1.0;
//			if constexpr (PROPERTIES == slabProfiles::ProfileProperties::MULTIPLYR) factor = radiusSamplesInternal_[i];
//			ProfileValues<slabProfiles::ProfileType::R>()[i] = factor * static_cast<Real>(rValues[i]);
//			ProfileValues<slabProfiles::ProfileType::T>()[i] = factor * static_cast<Real>(tValues[i]);
//		}
//
//	}
//
//	template class MCMLProfile<slabProfiles::ProfileProperties::NONE>;
//	template class MCMLProfile<slabProfiles::ProfileProperties::MULTIPLYR>;
//}
