#pragma once


struct dviwedi
{
	slabProfiles::Real get_highv0(double const alpha) const
	{
		const auto alphaoneminus = 1 - alpha;
		const auto potalphaoneminus = alphaoneminus * alphaoneminus;
		const auto tripalphaoneminus = potalphaoneminus * alphaoneminus;
		const auto highv0 = sqrt(3 * alphaoneminus) * (1 - (2 / 5)*alphaoneminus - 12 / 175 * potalphaoneminus - (2 / 125) * tripalphaoneminus - (166 / 67375) * tripalphaoneminus * alphaoneminus);

		return 1 / highv0;
	}

	slabProfiles::Real get_lowv0(double const alpha) const
	{

		slabProfiles::Real const pot_aplha = alpha * alpha;
		slabProfiles::Real const  trip_alpha = pot_aplha * alpha;

		slabProfiles::Real lowv0 = 1 - 2 * exp(-2 / alpha)*(1 + ((4 - alpha) / alpha)*exp(-2 / alpha) + ((24 - 12 * alpha + pot_aplha) / pot_aplha)*exp(-4 / alpha) +
			((512 - 384 * alpha + 72 * pot_aplha - 3 * trip_alpha) / trip_alpha)*exp(-6 / alpha));

		if (lowv0 == 1.0) lowv0 = 0.99999999;


		return static_cast<slabProfiles::Real>(1 / lowv0);
	}

	std::size_t static const amount_abs_scatter = 5;

	double abs_scatter_values[amount_abs_scatter] = {
		0.01,
		0.1,
		1.0,
		10.0,
		100.0
	};

	std::size_t static const amount_anisotropy = 3;

	double anisotropy[amount_anisotropy] = {
		0.0,
		0.5,
		0.9
	};


	std::size_t static const amount_photons = 6;
	int photons[amount_photons]{
		500,
		2000,
		5000,
		10000,
		50000,
		100000,
	};


};

