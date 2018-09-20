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
	void update_direction(photonstruct* photon, material const* mat)  override;
	double cal_stepsize(photonstruct* photon, material const* mat)  override;
	void run(const std::string mcml_path) override;
	double directional_distribution_dwi(double const v0, double const wz);
	double directional_distribution_hg(double const g, double const wz);
	void cal_absorption(photonstruct * photon, material const * mat_) const override;
	bool is_hit(photonstruct * photon, material const * mat_) override;
};

