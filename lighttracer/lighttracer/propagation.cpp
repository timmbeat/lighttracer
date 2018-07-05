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

double propagation::sign(double const x)
{
	if (x == 0) return 0;
	return x > 0 ? 1 : -1;
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
		//TODO: THIS LITTLE ONE HERE WAS THE MISTAKE ALL ALONG
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

void propagation::write_to_logfile(output * out_, material const *  mat_, std::string const path, std::string const csv_path) const
{
	auto mc = mcml::MCMLParser("C://Users//Tim//Documents//WISE17_18//Bachelor//FirstPrototype//mcml_test//bli.txt");
	std::ofstream fout(path, std::ofstream::trunc);
	std::ofstream ccout(csv_path, std::ofstream::trunc);


	auto bla = mc.get_ra();
	std::stringstream csvout;
	csvout << std::setw(11) << std::left << "weight_me" << std::setw(11) << std::left << "weight_mcml" << " ir";
	ccout << csvout.str() << std::endl;
	csvout.str("");
	std::stringstream output;
	auto counter = 0.0;
	///////////////////
	//Scaling/////////
	//////////////////
	auto scale1 = 2.0*slabProfiles::pi<double>()*out_->delr*out_->delr*mat_->num_photons;
	

	for (auto j = 0; j < out_->bins_r; j++)
	{
		for (auto i = 0; i < out_->bins_a; i++)
		{
			output.str("");
			csvout.str("");
			
				csvout << std::setw(10) << std::left << out_->Rd_r[j][i] * 1.0 / ((j + 0.5)*scale1) << " " << std::setw(11) << std::left << bla[j] << " " << j << std::endl;
				output << std::setw(10) << std::left << out_->Rd_r[j][i] * 1.0 / ((j + 0.5)*scale1) << " ir " << j << " ia " << i;
				ccout << csvout.str();
				fout << output.str();
				fout << std::endl;
			
			
			
					counter += out_->Rd_r[j][i];
				
			
			
		}
	}
	


	auto scale3 = 1.0 / static_cast<double>(mat_->num_photons);
	output.str("");
	output << "WEIGHT OF TRANSMITTANCE " << counter * scale3;
	fout << output.str();
	fout << std::endl;


}

//int main()
//{
//
//
//	auto mc = mcml::MCMLParser("C://Users//Tim//Documents//WISE17_18//Bachelor//FirstPrototype//mcml_test//bli.txt");
//
//	auto const photons = mc.get_numphotons();
//	auto const layer_mcml = mc.GetLayers();
//	auto const layer_0 = layer_mcml[0];
//	auto lay = layer(layer_0.eta_, layer_0.mua_, layer_0.mus_, layer_0.g_, propagation::getv0(layer_0.mus_ / (layer_0.mua_ + layer_0.mus_)));
//	auto material1 = material(photons, render.wth, 0.1, &lay, 0.999999999, false);
//	output out(mc.get_numr(), mc.get_numa(), mc.get_dr_(), 1);
//	std::cout << lay.reflec << std::endl;
//	propagation prop(material1, out);
//	for(auto i = 0; i < photons; i++)
//	{
//		photonstruct photon(1 - lay.reflec, false);
//		while(!photon.dead)
//		{
//			prop.trace(photon);
//			if (photon.weight < prop.get_material().wth && !photon.dead)
//			{
//				prop.roulette(photon);
//			}
//		}
//	}
//	prop.write_to_logfile();
//
//
//}


