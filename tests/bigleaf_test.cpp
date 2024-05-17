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

	string pfile = "tests/params/p_test_v2.ini";
	if (argc > 1) pfile = argv[1];

	int err = 0;

	pfate::BigLeafPatch sim(pfile);
	sim.init();
	sim.set_fapar(0.92);
	sim.simulate(2000, 2020);

	return err;

}


