#pragma once
#include "material.h"
#include "output.h"
#include "Photon.h"
#include <functional>
#include <random>
#include "mcml_parser.h"
#include "glm/glm.hpp"
static std::random_device rd;
static std::mt19937 gen(rd());
static std::uniform_real_distribution<double> dis(0.0, 1.0);

class propagation
{
	public:
	virtual ~propagation() = default;
	void virtual cal_absorption(photonstruct * photon, material const * mat) const;
	void move(photonstruct * photon);
	void trace(photonstruct * photon, output * out, material  const * mat);
	void update_arr_bucket(photonstruct const * photon, output * out);
	void roulette(photonstruct * photon, material const * mat) const;
	void virtual update_direction(photonstruct * photon, material const * mat);
	void virtual run(std::string mcml_path) = 0;

	double fresnel(double const uz, material const * mat) const;
	double sample_classical_distribution(double const mu_t);
	double virtual cal_stepsize(photonstruct * photon, material const * mat) = 0;;
	double classical_path_distribution(double mu_t, double stepsize) const;
	double get_threshould() const;
	
	
	glm::dvec2 calculate_scattering(double anisotropy);
	bool is_hit(photonstruct * photon, material const * mat);
	int ray_triangle(glm::dvec3 &point, glm::dvec3 &vector, std::vector<glm::dvec3> & plane);


	private:
	double const threshould_weight_ = 0.001;

};

static double random()
{
	return dis(gen);
}
