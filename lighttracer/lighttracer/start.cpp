#include "stdafx.h"
#include "dwivedi_sampling.h"

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
	dwivedi.run(12, 110, 0.99, 200000);
}
