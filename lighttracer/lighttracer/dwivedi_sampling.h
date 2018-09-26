#pragma once
#include "propagation.h"

class dwivedi_sampling : public propagation
{


	public:
	
	private:
	static double getv0(double const alpha);
	
	
	public:
	dwivedi_sampling() = default;
	~dwivedi_sampling() = default;
	void update_direction(photonstruct* photon, material const* mat) override;
	void run(const std::string mcml_path) override;
	void run(double mua, double mus, double anisotropy, std::size_t photons, std::size_t numr = 100, double delr = 0.001);
	double directional_distribution_dwi(double const v0, double const wz) const;
	double directional_distribution_hg(double const g, double const wz) const;
	void cal_absorption(photonstruct * photon, material const * mat_) const override;
	double dwivedi_path_distribution(double wz, double v0, double mu_t, double t) const;
	double cal_stepsize(photonstruct* photon, material const * mat) override;

	//bool is_hit(photonstruct * photon, material const * mat_) override;

};

