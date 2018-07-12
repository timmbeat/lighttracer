#pragma once
#include "material.h"
#include "Input.h"
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

	void virtual update_direction(photonstruct * photon, material const * mat) = 0;
	double virtual cal_stepsize(photonstruct * photon, material const * mat) = 0;
	void virtual run(const std::string mcml_path) = 0;


	void cal_absorption(photonstruct * photon, material const * mat_) const;
	void move(photonstruct * photon);
	
	void trace(photonstruct * photon, output * out, material  const * mat_);

	void update_arr_bucket(photonstruct const * photon, output * out_);
	void write_to_logfile(output * out_, material const *  mat_, std::string const path, std::string const csv_path, mcml::MCMLParser &mc) const;
	void roulette(photonstruct * photon, material const * mat_) const;
	double fresnel(double const uz, material const * mat_) const ;
	
	bool is_hit(photonstruct * photon, material const * mat_);


	int RayTriangle(glm::dvec3 &Point, glm::dvec3 &Vector, std::vector<glm::dvec3> & Plane);
	


};

static double random()
{
	return dis(gen);
}
