#pragma once
#include "propagation.h"

class dwivedi_sampling : propagation
{


	public:
	dwivedi_sampling() = default;
	~dwivedi_sampling() = default;
	private:
	static double getv0(double const alpha);
	public:
	void update_direction(photonstruct* photon, material const* mat)  override;
	double cal_stepsize(photonstruct* photon, material const* mat)  override;
	void run(const std::string mcml_path) override;
	double scattering_function(double const v0, double const wz);
	double scattering_function_hg(double const g, double const wz);
	void cal_absorption(photonstruct * photon, material const * mat_) const override;
	bool is_hit(photonstruct * photon, material const * mat_) override;
};

