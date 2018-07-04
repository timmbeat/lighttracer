#pragma once
#include <string>


#include "propagation.h"

class classic_sampling : propagation
{
	public:
	classic_sampling();
	~classic_sampling();

	private:
	static glm::dvec2 calculate_scattering(double anisotropy);
public:
	void run(const std::string mcml_path) override;
	void update_direction(photonstruct* photon, material const* mat) override;
	double cal_stepsize(photonstruct* photon, material const* mat) override;
};

