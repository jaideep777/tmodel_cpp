#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include <numeric>
#include <functional>
#include <algorithm>

#include "plantfate_patch.h"

using namespace std;

int main(int argc, char ** argv){

	if (argc < 4){
		cout << "Plant-FATE error: pf_splash needs three arguments. Syntax: ./pfate <params_file.ini> <start_year> <end_year>\n";
		return 1;
	}

	string pfile = argv[1];
	double y0 = stod(argv[2]), yf = stod(argv[3]);

	if (y0 > yf){
		cout << "Plant-FATE error: start year should be <= end year\n";
	}

	pfate::Patch sim(pfile);
	sim.climate_stream.use_precip_data = true;
	sim.init(y0, yf);
	sim.spinup();
	sim.simulate_coupled();
	sim.close();

	return 0;
}

