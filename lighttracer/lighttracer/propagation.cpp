#include "stdafx.h"
#include "propagation.h"
#include "config.h"
#include <fstream>
#include <random>
#include <iostream>
#include "mcml_parser.h"
std::random_device rd;
std::mt19937 gen(rd());
std::uniform_real_distribution<double> dis(0.0, 1.0);
const double m_pi = 3.14159265358979323846;
std::ofstream fout("logfile.txt", std::ofstream::app);
std::ofstream ccout("output.csv", std::ofstream::app);

extern struct dviwedi dvi;
extern struct renderoptions render;

propagation::~propagation()
= default;
//triple checkt
void propagation::cal_absorption(photonstruct& photon) const
{
	auto const tmp = photon.weight * mat_.matproperties->absorption / (mat_.matproperties->absorption + mat_.matproperties->scattering);
	photon.weight -= tmp;
}

//triple checkt...
double propagation::cal_stepsize(photonstruct& photon) const
{
	double ran;
	do
	{
		ran = dis(gen);

	}
	while (ran <= 0.0);
	auto const ret = -log(dis(gen)) / (mat_.matproperties->absorption + mat_.matproperties->scattering);

	return ret;
}

double propagation::sign(double const x)
{
	if (x == 0) return 0;
	return x > 0 ? 1 : -1;
}

void propagation::move(photonstruct& photon)
{
	photon.position.x = photon.position.x + photon.direction.x*photon.step;
	photon.position.y = photon.position.y + photon.direction.y*photon.step;
	photon.position.z = photon.position.z + photon.direction.z*photon.step;
}

//Triple checkt..
glm::dvec2 propagation::calculate_scattering() const
{
	double scattering;
	if (mat_.matproperties->anisotropy == 0.0)
	{

		scattering = 2 * dis(gen) - 1;
	}
	else
	{
		auto const g = mat_.matproperties->anisotropy;
		auto const temp = (1 - g * g) / (1 - g + 2 * g * dis(gen));
		scattering = (1 + g * g - temp * temp) / (2 * g);
	}

	auto const azimuthal = 2 * m_pi * dis(gen);

	return glm::dvec2(scattering, azimuthal);
}

//triple checkt..
void propagation::update_direction(struct photonstruct& photon) const
{
	auto const angles = calculate_scattering();
	auto const sint = sqrt(1 - angles.x * angles.x);
	auto const cosp = cos(angles.y);
	double sinp;
	auto const ux = photon.direction.x;
	auto const uy = photon.direction.y;
	auto const uz = photon.direction.z;
	if (angles.y < m_pi)
	{
		sinp = sqrt(1.0 - cosp * cosp);
	}
	else
	{
		sinp = -sqrt(1.0 - cosp * cosp);
	}

	if (fabs(photon.direction.z) > mat_.angle_threshold)
	{
		photon.direction.x = sint * cosp;
		photon.direction.y = sint * sinp;
		photon.direction.z = sign(uz)*angles.x;
	}
	else
	{
		auto const temp = sqrt(1.0 - uz * uz);
		photon.direction.x = sint * (ux*uz*cosp - uy * sinp) / temp + ux * angles.x;

		photon.direction.y = sint * (uy*uz*cosp + ux * sinp) / temp + uy * angles.x;

		photon.direction.z = -sint * cosp*temp + uz * angles.x;
	}


}

double propagation::dwivedi() const
{
	auto const y = dis(gen);
	auto const x = (mat_.matproperties->v0 - 1) / (mat_.matproperties->v0 + 1);
	auto const wz_new = mat_.matproperties->v0*pow((mat_.matproperties->v0 + 1), y)*(pow((mat_.matproperties->v0 - 1), y) + pow((mat_.matproperties->v0 + 1), y));
	auto const tmp = pow(x, y);
	auto const wz = mat_.matproperties->v0 - (mat_.matproperties->v0 + 1)*tmp;
	return (-log(1 - dis(gen))) / ((1 - wz / mat_.matproperties->v0)*(mat_.matproperties->scattering + mat_.matproperties->absorption));
}

double propagation::getv0(double const alpha)
{
	return alpha > 0.56 ? dvi.get_highv0(alpha) : dvi.get_lowv0(alpha);
}


void propagation::roulette(photonstruct &photon) const
{
	if(photon.weight == 0.0)
	{
		photon.dead = true;
	}
	else if (dis(gen) < mat_.rusroul)
	{
		photon.weight /= mat_.rusroul;
	}
	else
	{
		photon.dead = true;
	}
}



bool propagation::is_hit(photonstruct& photon)
{
	//TODO: Test for both sampling sampling_scheme
	
	if(photon.sleft == 0.0)
	{
		if (mat_.dwivedi)
		{
			photon.step = dwivedi();
		}
		else
		{
			photon.step = cal_stepsize(photon);
		}
	}
	else
	{
		photon.step = photon.sleft;
		photon.sleft = 0.0;
	}


	if (photon.direction.z < 0.0)
	{
		auto const s1 = photon.position.z / photon.direction.z;
		//No check for lower boundery, because there is none
		if (s1 < photon.step && photon.direction.z != 0.0)
		{
			photon.sleft = (photon.step - s1)*(mat_.matproperties->absorption + mat_.matproperties->scattering);
			photon.step = s1;

			return true;
		}
	}

	return false;
}


void propagation::trace(photonstruct &photon)
{
	//Specular reflectance here, but not necassary because that would only decrement
	//the weight by constant


	if (is_hit(photon))
	{
		move(photon);
		auto const cos90D = 1 / 1000000;
		auto const uz = -photon.direction.z;
		double r;
		if (uz > mat_.angle_threshold)
		{
			r = (1 - mat_.matproperties->refrac / 1 + mat_.matproperties->refrac);
			r *= r;
		}
		else if (uz < cos90D)
		{
			r = 1.0;
		}
		else
		{
			double sa1, sa2;
			double ca2;

			sa1 = sqrt(1 - uz * uz);
			sa2 = uz * sa1;

			if (sa2 >= 1.0)
			{
				r = 1.0;
			}
			else
			{
				double cap, cam;
				double sap, sam;
				ca2 = sqrt(1 - sa2 * sa2);

				cap = uz * ca2 - sa1 * sa2;
				cam = uz * ca2 + sa1 * sa2;
				sap = sa1 * ca2 + uz * sa2;

				sam = sa1 * ca2 - uz * sa2;
				r = 0.5 * sam*sam*(cam*cam + cap * cap) / (sap*sap*cam*cam);
			}

		}
		/*const auto angleofincidence = acos(abs(photon.direction.z));
		const auto angleoftransmission = asin(mat_.matproperties->refrac*sin(angleofincidence));
		auto const a = sin(angleofincidence - angleoftransmission);
		auto const b = sin(angleofincidence + angleoftransmission);

		auto const smtmp = a * a;
		auto const sptmp = b * b;
		auto const tmtmp = a / sqrt(1 - a * a);
		auto const tptmp = b / sqrt(1 - b * b);

		const auto internalreflactance = 1 / 2 * smtmp / sptmp +
			(tmtmp*tmtmp) / (tptmp*tptmp);*/

		if (dis(gen) <= r)
		{
			photon.direction.z *= -1;
			resume(photon);

		}
		else
		{
			update_arr_bucket(photon);
			//If photon is out of the tissue
			photon.dead = true;
			//updateBuckets(photon, res.trans);
		}
	}
	else
	{
		resume(photon);
	}

}

void propagation::resume(photonstruct & photon)
{
	move(photon);
	cal_absorption(photon);
	update_direction(photon);
	if (photon.weight < mat_.wth && !photon.dead)
	{
		roulette(photon);
	}
	else if (photon.weight < 0.0)
	{
		photon.dead = true;
	}
	if (!photon.dead)
	{
		trace(photon);
	}
}


void propagation::update_arr_bucket(photonstruct &photon)
{
	auto const r = sqrt(photon.position.x * photon.position.x + photon.position.y * photon.position.y);
	auto ir = static_cast<int>(r / out_.delr);
	if (ir > out_.ARRSIZE - 1) ir = out_.ARRSIZE - 1;

	//auto ia = static_cast<int>(acos(-photon.direction.z) / (2 * m_pi*(ir + 0.5)*(out_.delr * out_.delr)));
	auto ia = static_cast<int>(acos(-photon.direction.z) / out_.dela);

	if (ia > out_.ARRSIZE - 1) ia = out_.ARRSIZE - 1;
	out_.ref[ir][ia] += photon.weight*(1 - mat_.matproperties->reflec);

}

void propagation::write_to_logfile() const
{
	auto mc = mcml::MCMLParser("C://Users//Tim//Documents//WISE17_18//Bachelor//FirstPrototype//mcml_test//bli.txt");

	std::string sampling_scheme = "";
	if (mat_.dwivedi)
	{
		sampling_scheme = " DWIVDI ";
	}
	else
	{
		sampling_scheme = " STANDART ";
	}

	std::cout << mc.get_numr() << std::endl;
	std::cout << mc.get_numa() << std::endl;
	auto bla = mc.get_ra();
	std::stringstream csvout;
	csvout << "weight_mcml " << "weight_me " << "ir";
	ccout << csvout.str() << std::endl;
	csvout.str("");
	std::stringstream output;
	auto counter = 0.0;
	for (auto j = 0; j < out_.ARRSIZE; j++)
	{
		for (auto i = 0; i < out_.ARRSIZE; i++)
		{
			output.str("");
			csvout.str("");
			if (out_.ref[j][i] != 0.0)
			{
				csvout << out_.ref[j][i] << " " << bla[j][0] << " " << j << std::endl;
				output << out_.ref[j][i] << " ir " << j << " ia " << i;
				ccout << csvout.str();
				fout << output.str();
				fout << std::endl;
			}

			counter += out_.ref[j][i];
		}
	}
	


	auto scale2 = static_cast<double>(1.0 / mat_.num_photons);

	output << "WEIGHT OF TRANSMITTANCE " << counter * scale2 << sampling_scheme;
	fout << output.str();
	fout << std::endl;


}

material propagation::get_material() const
{
	return mat_;
}