#include "stdafx.h"
#include "propagation.h"
#include "mcml_parser.h"
#include "constants.h"
#include <functional>
#include "output.h"

void propagation::update_direction(photonstruct* photon, material const* mat)
{
	auto const angles = calculate_scattering(mat->matproperties->anisotropy);
	//auto const angles = calculate_scattering();
	auto const sint = sqrt(1 - angles.x * angles.x);
	auto const cosp = cos(angles.y);
	double sinp;
	auto const ux = photon->direction.x;
	auto const uy = photon->direction.y;
	auto const uz = photon->direction.z;
	if (angles.y < slabProfiles::pi<double>())
	{
		sinp = sqrt(1.0 - cosp * cosp);
	}
	else
	{
		sinp = -sqrt(1.0 - cosp * cosp);
	}

	if (fabs(photon->direction.z) > slabProfiles::cos_1<double>())
	{
		photon->direction.x = sint * cosp;
		photon->direction.y = sint * sinp;
		photon->direction.z = glm::sign(uz)*angles.x;
	}
	else
	{
		auto const temp = sqrt(1.0 - uz * uz);
		photon->direction.x = sint * (ux*uz*cosp - uy * sinp) / temp + ux * angles.x;

		photon->direction.y = sint * (uy*uz*cosp + ux * sinp) / temp + uy * angles.x;

		photon->direction.z = -sint * cosp*temp + uz * angles.x;
	}
}

glm::dvec2 propagation::calculate_scattering(double const anisotropy)
{
	auto const g = anisotropy;
	double scattering;
	if (g == 0.0)
	{

		scattering = 2 * random() - 1;
	}
	else
	{
		auto const temp = (1 - g * g) / (1 - g + 2 * g * random());
		scattering = (1 + g * g - temp * temp) / (2 * g);
	}

	auto const azimuthal = 2 * slabProfiles::pi<double>() * random();

	return glm::dvec2(scattering, azimuthal);
}

double propagation::sampleClassicalDistribution(double const mu_t)
{
	double ran;
	do
	{
		ran = random();

	}
	while (ran <= 0.0);
	auto const ret = -log(ran) / (mu_t);

	return ret;
}

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

double propagation::classical_path_distribution(double mu_t, double stepsize)
{
	return mu_t * exp(-mu_t * stepsize);
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
/*
 * Maybe Remove this. 
 * This was intended for Camera and Light implementations
 */
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

double propagation::getThreshould() const
{
	return threshould_weight_;
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
	else if(uz < slabProfiles::cos90_d<double>())
	{
		r = 1.0;
	}
	else
	{
		auto const sa1 = sqrt(1 - uz * uz);
		auto const sa2 = mat_->matproperties->refrac * sa1;
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

/*
 * Write into the Buckets 
 * 
 */
void propagation::update_arr_bucket(photonstruct const * photon, output * out_)
{
	auto ir = static_cast<int>(sqrt(photon->position.x * photon->position.x + photon->position.y * photon->position.y)
							   / out_->delr);
	if (ir > out_->bins_r - 1) ir = out_->bins_r - 1;

	auto ia = static_cast<int>((acos(-photon->direction.z) / out_->dela));
	if (ia > out_->bins_a - 1)
	{
		ia = out_->bins_a - 1;

	}
	out_->Rd_r[ir][ia] += photon->weight;

}
