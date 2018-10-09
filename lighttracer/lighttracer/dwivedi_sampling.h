#pragma once
#include "propagation.h"

class dwivedi_sampling : public propagation
{


	public:



	dwivedi_sampling() = default;
	~dwivedi_sampling() = default;

	/* Method to start Dwivedi run by using a mcml output file. This one is much faster than the other run Method
	 * @parameter mcml_path : path to the mcml output file
	 *
	 */
	void run(const std::string mcml_path) override;


	/* Method to start own Dwivedi run without using a mcml output file. This will actually conduct 2 runs. One 
	 * for classical Sampling and the Dwivedi one. 
	 * @parameter mua : absorption coefficent
	 * @parameter mis : scattering coefficient
	 * @parameter anisotropy : anisotropy
	 * @parameter photons : amount of photons used
	 * @parameter numr : amount of buckets [can be empty]
	 * @parameter delr : size of buckets [can be empty]
	 * @return : returns nothing
	 */
	void run(double mua, double mus, double anisotropy, std::size_t photons, std::size_t numr = 100, double delr = 0.001);
	void data_run();


	double cal_stepsize(photonstruct* photon, material const* mat) override;

	private:


	static double getv0(double const alpha);
	void update_direction(photonstruct* photon, material const* mat) override;
	double directional_distribution_dwi(double const v0, double const wz) const;
	double directional_distribution_hg(double const g, double const theta) const;
	void cal_absorption(photonstruct * photon, material const * mat) const override;
	double dwivedi_path_distribution(double wz, double v0, double mu_t, double t) const;


	///////DATA EXTRACTION

	
};

