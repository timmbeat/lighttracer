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
};

