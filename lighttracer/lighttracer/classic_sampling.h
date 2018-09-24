#pragma once
#include <string>


#include "propagation.h"

class classic_sampling : public propagation
{
	public:
	classic_sampling();
	~classic_sampling();

	private:
	static glm::dvec2 calculate_scattering(double anisotropy);
public:
	void run(const std::string mcml_path) override;
	double cal_stepsize(photonstruct* photon, material const* mat) override;
};

