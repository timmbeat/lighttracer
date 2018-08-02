#include "stdafx.h"
#include "propagation.h"
#include "config.h"
#include <fstream>
#include <random>
#include <iostream>
#include "mcml_parser.h"
#include "constants.h"
#include <functional>
#include <iomanip>

struct renderoptions render;
//triple checkt

void propagation::cal_absorption(photonstruct * photon, material const * mat_) const
{
	auto const tmp = photon->weight * mat_->matproperties->absorption / (mat_->matproperties->absorption + mat_->matproperties->scattering);
	
	photon->weight -= tmp;
}


void propagation::move(photonstruct * photon)
{
	photon->position.x = photon->position.x + photon->direction.x*photon->step;
	photon->position.y = photon->position.y + photon->direction.y*photon->step;
	photon->position.z = photon->position.z + photon->direction.z*photon->step;
}

void propagation::roulette(photonstruct  * photon, material const * mat_) const
{
	if(photon->weight == 0.0)
	{
		photon->dead = true;
	}
	else if (random() < mat_->rusroul)
	{
		photon->weight /= mat_->rusroul;
	}
	else
	{
		photon->dead = true;
	}
}



bool propagation::is_hit(photonstruct * photon, material const * mat_)
{
	//TODO: Test for both sampling sampling_scheme
	
	if(photon->sleft == 0.0)
	{
		photon->step = cal_stepsize(photon, mat_);
	}
	else
	{
		photon->step = photon->sleft / (mat_->matproperties->absorption + mat_->matproperties->scattering);
		photon->sleft = 0.0;
	}


	if (photon->direction.z < 0.0)
	{
		auto const s1 = -photon->position.z / photon->direction.z;
		//No check for lower boundery, because there is none
		if (s1 < photon->step && photon->direction.z != 0.0)
		{
			photon->sleft = (photon->step - s1)*(mat_->matproperties->absorption + mat_->matproperties->scattering);
			photon->step = s1;

			return true;
		}
	}

	return false;
}

void propagation::trace(photonstruct * photon, output * out, material  const * mat_)
{

	if (is_hit(photon, mat_))
	{

		move(photon);
		auto const uz = photon->direction.z;
		double r;

		if(-uz < slabProfiles::cos90_d<mcml::Real>())
		{
			r = 1.0;
		}
		else r = fresnel(-uz, mat_);

		if (random() <= r)
		{
			photon->direction.z = -uz;
		}
		else
		{
			update_arr_bucket(photon, out);
			//If photon is out of the tissue
			/*auto Tri = std::vector<glm::dvec3>{ glm::dvec3(-1.0, -1.0, -2.0) , glm::dvec3(1.0,-1.0, -2.0), glm::dvec3(0.0, 1.0, -2.0) };
			if(RayTriangle(photon->position, photon->direction, Tri) == 1)
			{
				
				out->kam_counter++;
				
				if( out->kam_counter % 500 == 0)
				{
					std::cout << out->kam_counter;
					std::cin.get();
				}
			}
*/
			photon->dead = true;
		}
	}
	else
	{
		move(photon);
		cal_absorption(photon, mat_);
		update_direction(photon, mat_);
	}

}

int propagation::RayTriangle(glm::dvec3 &Point, glm::dvec3 &Vector, std::vector<glm::dvec3> &Plane)
{
	auto const a = Plane[0];
	auto const b = Plane[1];
	auto const c = Plane[2];
		  
	auto const xamxb = a.x - b.x; //---a---
	auto const xamxc = a.x - c.x; //---d---
	auto const xamxe = a.x - Point.x; //---j---
		 
	auto const yamye = a.y - Point.y; //--k--
	auto const yamyc = a.y - c.y; //--e--
	auto const yamyb = a.y - b.y; //--b--
		 
	auto const zamze = a.z - Point.z; //--l--
	auto const zamzc = a.z - c.z; //--f--
	auto const zamzb = a.z - b.z; //--c--

	auto const M = xamxb * (yamyc*Vector.z - Vector.y*zamzc) + yamyb * (Vector.x*zamzc - xamxc * Vector.z) + zamzb * (xamxc*Vector.y - yamyc * Vector.x);


	auto t = - (zamzc *(xamxb*yamye - xamxe * yamyb) + yamyc * (xamxe*zamzb - xamxb * zamze) + xamxc * (yamyb*zamze - yamye * zamzb)) / M;

	if (t < 0) return 0;
	
	auto gamma = (Vector.z *(xamxb*yamye - xamxe * yamyb) + Vector.y*(xamxe*zamzb - xamxb * zamze) + Vector.x*(yamyb*zamze - yamye * zamzb)) / M;

	if (gamma < 0 || gamma > 1) return false;

	auto beta = (xamxe*(yamyc*Vector.z - Vector.y*zamzc) + yamye * (Vector.x*zamzc - xamxc * Vector.z) + zamze * (xamxc*Vector.y - yamyc * Vector.x)) / M;

	if (beta < 0 || beta > 1 - gamma) return false;

	return true;
	
}

double propagation::fresnel(double const uz, material const * mat_) const 
{
	//Because Vaccum has a refractive index of 1.0
	double r;
	if(mat_->matproperties->refrac == 1.0)
	{
		r = 0.0;
	}
	else if(uz > slabProfiles::cos_1<double>())
	{
		r = (mat_->matproperties->refrac-1) / (mat_->matproperties->refrac+1);
		r *= r;
	}
	//Insert Constant here
	else if(uz < slabProfiles::cos90_d<double>())
	{
		r = 1.0;
	}
	else
	{
		auto const sa1 = sqrt(1 - uz * uz);
		auto const sa2 = mat_->matproperties->refrac * sa1;
		/*
		 * OLD ONE: auto const sa2 = uz * sa1..
		 */
		if (sa2 >= 1.0)
		{
			r = 1.0;
		}
		else
		{

			auto const ca2 = sqrt(1 - sa2 * sa2);

			auto const cap = uz * ca2 - sa1 * sa2;
			auto const cam = uz * ca2 + sa1 * sa2;
			auto const sap = sa1 * ca2 + uz * sa2;

			auto const sam = sa1 * ca2 - uz * sa2;
			r = 0.5 * sam*sam*(cam*cam + cap * cap) / (sap*sap*cam*cam);
		}

	}
	
	return r;
}

void propagation::update_arr_bucket(photonstruct const * photon, output * out_)
{
	auto ir = static_cast<int>(sqrt(photon->position.x * photon->position.x + photon->position.y * photon->position.y)
							   / out_->delr);
	if (ir > out_->bins_r - 1) ir = out_->bins_r - 1;

	//auto ia = static_cast<int>(acos(-photon.direction.z) / (2 * m_pi*(ir + 0.5)*(out_.delr * out_.delr)));
	auto ia = static_cast<int>((acos(-photon->direction.z) / out_->dela));
	if (ia > out_->bins_a - 1)
	{
		ia = out_->bins_a - 1;

	}
	out_->Rd_r[ir][ia] += photon->weight;

}
