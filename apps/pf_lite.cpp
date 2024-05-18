#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include <numeric>
#include <functional>
#include <algorithm>

#include "bigleaf_patch.h"

using namespace std;

int is_equal(const vector<double>& v1, const vector<double>& v2, double tol=1e-6){
	bool b = true;
	for (int i=0; i<v1.size(); ++i){
		cout << "Comparing v["<<i<<"]: " << v1[i] << " " << v2[i] << '\n';
		b = b & (fabs(v1[i]-v2[i])< tol);
	}
	return b? 0:1;
}

int main(int argc, char ** argv){

	if (argc < 4){
		cout << "Plant-FATE Lite error: pf_lite needs three arguments. Syntax: ./pf_lite <params_file.ini> <start_year> <end_year>\n";
		return 1;
	}

	string pfile = argv[1];
	double y0 = stod(argv[2]), yf = stod(argv[3]);

	if (y0 > yf){
		cout << "Plant-FATE Lite error: start year should be <= end year\n";
	}

	pfate::BigLeafPatch sim(pfile);
	sim.climate_stream.use_precip_data = true;
	sim.init_co2(414);
	sim.set_fapar(0.92);
	sim.init(2000, 2010.999);
	sim.simulate();

	return 0;
}

