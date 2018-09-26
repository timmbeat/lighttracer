#include "stdafx.h"
#include "classic_sampling.h"
#include "dwivedi_sampling.h"
#include <iostream>

/**
 * Created by Tim Mend.
 * See folder "mcml_examples" for possible input files
 * 
 */



int main()
{

	/*classic_sampling classic{};
	classic.run("mcml_examples/example_2.txt");
*/
	dwivedi_sampling dwivedi{};
	//dwivedi.run("mcml_examples/example_10.txt");
	dwivedi.run(75, 83, 0.89, 100000);
}
