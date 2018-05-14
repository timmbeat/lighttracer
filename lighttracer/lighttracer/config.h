#pragma once


//Skin Layer
struct skin
{
	//Absorption coefficient for human skin according to https://www.spiedigitallibrary.org/journals/journal-of-biomedical-optics/volume-17/issue-9/090901/Optical-properties-of-human-skin/10.1117/1.JBO.17.9.090901.full
	const double absorption = 0.25;
	//Skin Scattering coefficient according to https://www.sciencedirect.com/science/article/pii/S0022202X15414836
	const double scattering = 0.78;
	const double anistropy = 0.6;
	//Human skin refraction index according to http://bmlaser.physics.ecu.edu/literature/2006%2003_human%20skin%20index.pdf
	const double refrac = 1.3686;
	const double alpha = 0.969;

}skin;

struct renderoptions
{
	//minimal weight which photons can have
	const double wth = 0.0001;
	//Amount of photons
	const int photons = 1000;
}renderoptions;

struct dviwedi
{
	
	double get_highv0(double const alpha) const
	{
		const auto alphaoneminus = 1 - alpha;
		const auto potalphaoneminus = alphaoneminus * alphaoneminus;
		const auto tripalphaoneminus = potalphaoneminus * alphaoneminus;
		const auto highv0 = sqrt(3 * alphaoneminus) * (1 - (2 / 5)*alphaoneminus - 12 / 175 * potalphaoneminus - (2 / 125) * tripalphaoneminus + (166 / 67375) * tripalphaoneminus * alphaoneminus);

		return 1 / highv0;
	}

	double get_lowv0(double const alpha) const
	{

		auto const pot_aplha = alpha * alpha;
		auto const  trip_alpha = pot_aplha * alpha;
		auto const lowv0 = 1 - 2 * exp(-2 / alpha)*(1 + ((4 - alpha) / alpha)*exp(-2 / alpha) + ((24 - 12 * alpha + pot_aplha) / pot_aplha)*exp(-4 / alpha) +
			((512 - 384 * alpha + 72 * pot_aplha - 3 * trip_alpha) / 3 * trip_alpha)*exp(-6 / alpha));

		return 1 / lowv0;
	}
}dvi;