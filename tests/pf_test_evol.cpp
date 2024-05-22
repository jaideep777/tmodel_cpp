#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include <numeric>
#include <functional>
#include <algorithm>

#include "plantfate_patch.h"

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

	string pfile = "tests/params/p_test_v2_evol2.ini";
	if (argc > 1) pfile = argv[1];

	int err = 0;

	for (int i=0; i<1; ++i){
		pfate::Patch sim(pfile);
		// sim.expt_dir = sim.expt_dir + "_414ppm";
		sim.E.init_co2(414);
		sim.init(2000, 2200);
		sim.simulate();

		vector<double> ba = sim.props.species.basal_area_vec;
		for (auto& b : ba) b*=1e4;
		cout << setprecision(10) << "Basal areas [m2/Ha]: " << ba << '\n';

		err = is_equal(ba, {13.90109726, 0, 15.36616065, 0});

		sim.close();
	}

	return err;


}


