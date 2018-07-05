#include "stdafx.h"
#include "classic_sampling.h"
#include "dwivedi_sampling.h"

int main()
{

	classic_sampling classic{};
	classic.run("mcml_out.txt");

	dwivedi_sampling dwivedi{};
	dwivedi.run("mcml_out.txt");

}
