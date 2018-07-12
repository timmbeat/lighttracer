#include "stdafx.h"
#include "classic_sampling.h"
#include "dwivedi_sampling.h"
#include <iostream>

/**
 * Created by Tim Mend.
 * You can use as input either
 * mcml_examples/example_1.txt or mcml_examples/example_2.txt
 */



int main()
{

	/*classic_sampling classic{};
	classic.run("mcml_examples/example_2.txt");
*/
	dwivedi_sampling dwivedi{};
	dwivedi.run("mcml_examples/example_2.txt");

}
