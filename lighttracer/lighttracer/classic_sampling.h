#pragma once
#include <string>


#include "propagation.h"

class classic_sampling : public propagation
{
	public:
	classic_sampling();
	~classic_sampling();

	private:
	public:
	void run(const std::string mcml_path) override;
	double cal_stepsize(photonstruct* photon, material const* mat) override;
};

