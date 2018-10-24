#include "stdafx.h"
#include "dwivedi_sampling.h"

/**
 * Created by Tim Mend.
 * See folder "mcml_examples" for possible input files
 * 
 */



int main()
{
	dwivedi_sampling dwivedi{};

	//Three different possibilites are aviable to start a dwivedi run.
	//First: use a mcml output file	
	dwivedi.run("mcml_examples/abs10sca0.1g0pho15000.txt");
	
	
	
	//Second: use your own valus
	//dwivedi.run(0.01, 100, 0.0, 5000);



	//Third: Specifiy some data in config.h and run it ===> This will take some time.
	//dwivedi.data_run();
}
