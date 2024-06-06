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

	string pfile = "tests/params/p_test_v2_gfguy.ini";
	if (argc > 1) pfile = argv[1];

	int err = 0;

	pfate::BigLeafPatch sim(pfile);
	sim.climate_stream.use_precip_data = true;
	sim.init_co2(414);
	sim.set_fapar(0.92);
	// sim.init(1995, 1995.3);
	// sim.init(1995, 2010.999);
	// sim.init(0, 109.5728);
	sim.init(0, 110);

	sim.spinup();
	cout << setprecision(12) << sim.soil_env.state.wn << "\n";
	if (fabs(sim.soil_env.state.wn - 773.199185007) > 1e-6) return 1; // uncoupled run

	sim.simulate();
	cout << setprecision(12) << sim.soil_env.state.wn << "\n";
	// if (fabs(sim.soil_env.state.wn - 740.283568986) > 1e-6) return 1; // uncoupled run
	// if (fabs(sim.soil_env.state.wn - 740.032327052) > 1e-6) return 1; // coupled run exchanging psi_m and trans, but still using file ppfd
	// if (fabs(sim.soil_env.state.wn - 740.40970978) > 1e-6) return 1; // coupled run exchanging psi_m and trans, ppfd set from splash
	if (fabs(sim.soil_env.state.wn - 740.409737135) > 1e-6) return 1; // coupled run exchanging psi_m and trans, ppfd set from splash, using days as time unit, note this unit change leads to a change in output in 4th decimal... FIXME: this could be a result of using avg year length earlier, but investigate

	return err;

}

