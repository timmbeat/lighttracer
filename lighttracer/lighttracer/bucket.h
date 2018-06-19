#pragma once


struct bucket
{
	
	double x;
	double y;
	double z;
	bucket(double const x, double const y, double const z, double const weight) : x(x), y(y), z(z), weight(weight), hits(1){ }
	double weight;
	int hits;
};

struct sort_bucket
{
	bool operator() (const bucket& buck, const bucket& buck2)const
	{
		return(buck.hits > buck2.hits);
	}
};