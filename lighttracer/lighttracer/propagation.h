#pragma once
#include "material.h"
#include "Input.h"
#include "Photon.h"




class propagation
{
	public:
	propagation(material const mat, output out) : mat_(mat), out_(out)
	{
	}
	~propagation();


	void cal_absorption(photonstruct &photon) const;
	void move(photonstruct &photon);
	void update_direction(photonstruct &photon)const;
	void trace(photonstruct &photon);
	void update_arr_bucket(photonstruct &photon);
	void write_to_logfile() const;
	void roulette(photonstruct &photon) const;
	double fresnel(double uz) const ;
	double cal_stepsize(photonstruct &photon) const;
	double dwivedi() const;

	static double getv0(double alpha);
	static double sign(double x);

	bool is_hit(photonstruct &photon);

	material get_material() const ;
	glm::dvec2 calculate_scattering() const;




	private:
	material mat_;
	output out_;

};


