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

	string pfile = "tests/params/p_test_v2_gfguy.ini";
	if (argc > 1) pfile = argv[1];

	int err = 0;

	pfate::Patch sim(pfile);
	sim.climate_stream.use_precip_data = true;
	// sim.init(1995.001, 1995.3);
	sim.init(0, (2050-1995+1)*365);

	sim.spinup();
	cout << setprecision(12) << sim.soil_env.state.wn << " " << sim.soil_env.dsoil.psi_m << "\n";
	if (fabs(sim.soil_env.state.wn - 773.199185007) > 1e-6) return 1; // uncoupled run

	sim.simulate_coupled();
	cout << setprecision(12) << sim.soil_env.state.wn << " " << sim.soil_env.dsoil.psi_m << "\n";

	if (fabs(sim.soil_env.state.wn - 739.975303616) > 1e-6) return 1; // coupled run exchanging psi_m and trans, ppfd set from splash, 7 day veg update

	sim.close();

	return err;

}

