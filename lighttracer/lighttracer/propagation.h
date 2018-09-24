#pragma once
#include "material.h"
#include "output.h"
#include "Photon.h"
#include <functional>
#include <random>
#include "mcml_parser.h"

static std::random_device rd;
static std::mt19937 gen(rd());
static std::uniform_real_distribution<double> dis(0.0, 1.0);

class propagation
{
	public:
	virtual ~propagation() = default;

	void virtual update_direction(photonstruct * photon, material const * mat);
	glm::dvec2 calculate_scattering(double anisotropy);
	double virtual cal_stepsize(photonstruct * photon, material const * mat) = 0;
	void virtual run(const std::string mcml_path) = 0;
	double sampleClassicalDistribution(double const mu_t);

	void virtual cal_absorption(photonstruct * photon, material const * mat_) const;
	void move(photonstruct * photon);
	
	void trace(photonstruct * photon, output * out, material  const * mat_);

	void update_arr_bucket(photonstruct const * photon, output * out_);
	void roulette(photonstruct * photon, material const * mat_) const;
	double fresnel(double const uz, material const * mat_) const ;

	virtual bool is_hit(photonstruct * photon, material const * mat_);

	double classical_path_distribution(double mu_t, double stepsize);
	int RayTriangle(glm::dvec3 &Point, glm::dvec3 &Vector, std::vector<glm::dvec3> & Plane);
	double getThreshould() const;


private:
double const threshould_weight_ = 0.001;

};

static double random()
{
	return dis(gen);
}
