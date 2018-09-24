#pragma once

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
};

