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
	//dwivedi.run("mcml_examples/Research/Skin.txt");
	//dwivedi.run(1.0, 10, 0.90, 10000000);
	dwivedi.data_run();
}
