/**
* @file   mcml_parser.h
* @author Sebastian Maisch <sebastian.maisch@uni-ulm.de>
* @date   2018.02.07
*
* @brief  Parser for MCML output files.
*/

#pragma once

//#include "depth_profile_generic.h"

#include <string>
#include <vector>
#include "constants.h"
#include <memory>

namespace mcml
{

	using Real = slabProfiles::Real;

	struct LayerInfo
	{
		Real eta_;
		Real mua_;
		Real mus_;
		Real g_;
		Real depth_;
	};

	//template<slabProfiles::ProfileProperties PROPERTIES> class MCMLProfile : public slabProfiles::DepthProfileGeneric
	//{
	//	public:
	//	MCMLProfile(Real alphap, Real sigmapt, Real g, Real eta, Real depth, const std::vector<Real>& radiusSamples,
	//				const std::vector<Real>& rValues, const std::vector<Real>& tValues);

	//	private:
	//	/** The internal radius samples used for this profile. */
	//	std::vector<Real> radiusSamplesInternal_;
	//};


	class MCMLParser
	{
		public:
		MCMLParser(const std::string& filename);

		const std::string& GetFilename() const
		{
			return filename_;
		}
		std::vector<Real> get_ra() const;
		const std::vector<LayerInfo>& GetLayers() const
		{
			return layers_;
		}
		Real GetEtaAbove() const
		{
			return etaTop_;
		}
		Real GetEtaBelow() const
		{
			return etaBelow_;
		}

		std::size_t get_numr() const;
		std::size_t get_numa() const;
		std::size_t get_numphotons() const;
		double get_dr_() const;
		double get_dz_() const;
		double get_diffuse_reflectance() const;
		// const std::vector<double>& GetRadiusValues() const { return rValues_; }
		// const std::vector<double>& GetRadialReflectionProfile() const { return rProfileR_; }
		// const std::vector<double>& GetRadialTransmissionProfile() const { return tProfileR_; }
		// const std::vector<double>& GetRadialReflectionProfileSmooth() const { return rProfileRSmooth_; }
		// const std::vector<double>& GetRadialTransmissionProfileSmooth() const { return tProfileRSmooth_; }

		/*template<slabProfiles::ProfileProperties PROPERTIES> std::unique_ptr<MCMLProfile<PROPERTIES>> GetProfile() const;
		template<slabProfiles::ProfileProperties PROPERTIES> std::unique_ptr<MCMLProfile<PROPERTIES>> GetProfileSmooth() const;*/
		void SmoothData();

		private:
		void ProcessInParams(std::ifstream& mcoIn);
		void ProcessRAT(std::ifstream& mcoIn);
		void ProcessFloatList(std::ifstream& mcoIn, std::vector<Real>& floats, Real factor = 1.0);
		void ProcessFloatList2(std::ifstream& mcoIn, std::vector<std::vector<Real>>& floats, Real factor = 1.0);
		void NextList2(std::size_t& i, std::size_t& j, const std::vector<std::vector<Real>>& floats);

		

		std::string SkipEmptyLines(std::ifstream& mcoIn);

		/** The filename of the mcml file. */
		std::string filename_;

		/** The number of photons used. */
		std::size_t numPhotons_;
		/** Size of a bin in z directions. */
		double dr_;
		/** Size of a radial bin. */
		double dz_;
		/** Number of bins in z directions. */
		std::size_t numBinsZ_;
		/** Number of radial bins. */
		std::size_t numBinsR_;
		/** number of angular bins. */
		std::size_t numBinsA_;

		/** Index of refraction on top. */
		Real etaTop_;
		/** Information on the layers. */
		std::vector<LayerInfo> layers_;
		/** Index of refraction below last layer. */
		Real etaBelow_;

		/** Specular reflectance. */
		double specularReflectance_;
		/** Diffuse reflectance. */
		double diffuseReflectance_;
		/** Absorbed fraction. */
		double absorbedFract_;
		/** Transmittance. */
		double transmittance_;

		/** Absorptions per layer. */
		std::vector<Real> absorptionsLayer_;
		/** Absorption by z. */
		std::vector<Real> absorptionsZ_;
		/** Reflectance by r. */
		std::vector<Real> rProfileR_;
		/** Reflectance by a. */
		std::vector<Real> rProfileA_;
		/** Transmittance by r. */
		std::vector<Real> tProfileR_;
		/** Transmittance by a. */
		std::vector<Real> tProfileA_;
		/** Absorption by r and z. */
		std::vector<std::vector<Real>> absorptionsRZ_;
		/** Reflectance by r and a. */
		std::vector<std::vector<Real>> rProfileRA_;
		/** Transmittance by r and a. */
		std::vector<std::vector<Real>> tProfileRA_;

		/** z values. */
		std::vector<Real> zValues_;
		/** r values. */
		std::vector<Real> rValues_;
		/** a values. */
		std::vector<Real> aValues_;

		/** r and t profiles (by r) smoothed */
		std::vector<Real> rProfileRSmooth_;
		std::vector<Real> tProfileRSmooth_;
	};


	/*template<slabProfiles::ProfileProperties PROPERTIES> std::unique_ptr<MCMLProfile<PROPERTIES>> MCMLParser::GetProfile() const
	{
		double sigmaps = layers_[0].mus_ * (1.0 - layers_[0].g_);
		double sigmapt = sigmaps + layers_[0].mua_;
		double alphap = sigmaps / sigmapt;
		return std::make_unique<MCMLProfile<PROPERTIES>>(alphap, sigmapt, layers_[0].g_, layers_[0].eta_, layers_[0].depth_, rValues_, rProfileR_, tProfileR_);
	}

	template<slabProfiles::ProfileProperties PROPERTIES> std::unique_ptr<MCMLProfile<PROPERTIES>> MCMLParser::GetProfileSmooth() const
	{
		if (rProfileRSmooth_.empty()) return GetProfile<PROPERTIES>();

		Real sigmaps = layers_[0].mus_ * (1.0 - layers_[0].g_);
		Real sigmapt = sigmaps + layers_[0].mua_;
		Real alphap = sigmaps / sigmapt;
		return std::make_unique<MCMLProfile<PROPERTIES>>(alphap, sigmapt, layers_[0].g_, layers_[0].eta_, layers_[0].depth_, rValues_, rProfileRSmooth_, tProfileRSmooth_);
	}
*/
	// namespace MCMLParserHelper {
	//     template<deepProfile::ProfileProperties PROPERTIES> std::unique_ptr<MCMLProfile<PROPERTIES>> GetProfile(const MCMLParser& parser)
	//     {
	//         double sigmaps = layers_[0].mus_ * (1.0 - layers_[0].g_);
	//         double sigmapt = sigmaps + layers_[0].mua_;
	//         double alphap = sigmaps / sigmapt;
	//         return std::make_unique<MCMLProfile<PROPERTIES>>(alphap, sigmapt, layers_[0].eta_, layers_[0].depth_, rValues_, rProfileR_, tProfileR_);
	//     }
	// 
	//     template<deepProfile::ProfileProperties PROPERTIES> std::unique_ptr<MCMLProfile<PROPERTIES>> GetProfileSmooth(const MCMLParser& parser)
	//     {
	//         if (rProfileRSmooth_.empty()) return GetProfile();
	// 
	//         double sigmaps = layers_[0].mus_ * (1.0 - layers_[0].g_);
	//         double sigmapt = sigmaps + layers_[0].mua_;
	//         double alphap = sigmaps / sigmapt;
	//         return std::make_unique<MCMLProfile<PROPERTIES>>(alphap, sigmapt, layers_[0].eta_, layers_[0].depth_, rValues_, rProfileRSmooth_, tProfileRSmooth_);
	//     }
	// }
}
