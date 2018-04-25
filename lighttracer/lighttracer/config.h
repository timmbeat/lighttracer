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
	double alpha = 0.959;
	double highv0 = sqrt(3 * (1 - alpha)) * (1 - 2 / 5*(1 - alpha) - 12 / 175 * pow((1 - alpha), 2) - 2 / 125 * pow((1 - alpha), 3) - 166 / 67375 * pow((1 - alpha), 4));
	double lowv0 = 1 - 2 * exp(-2 / alpha)*(1 + ((4 - alpha) / alpha)*exp(-2 / alpha) + ((24 - 12 * alpha + pow(alpha, 2)) / pow(alpha, 2))*exp(-4 / alpha) + 
		((512 - 384 * alpha + 72 * pow(alpha, 2) - 3 * pow(alpha, 3)) / pow(alpha, 3))*exp(-6 / alpha));
}dvi;