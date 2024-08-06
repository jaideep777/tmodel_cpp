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

	string pfile = "tests/params/p_test_v2.ini";
	if (argc > 1) pfile = argv[1];

	int err = 0;

	for (int i=0; i<1; ++i){
		pfate::Patch sim(pfile);

		// running with year as time unit
		double tpy = 1; // time units per year
		sim.config.time_unit = "years since 0000-01-00 0:0:0";

		// // running with day as time unit
		// double tpy = 365.2425;
		// sim.config.time_unit = "days since 0000-01-00 0:0:0";

		// translate all time invervals to new units
		sim.config.timestep = 0.04166666666666666666666*tpy;

		sim.E.init_co2(414);
		sim.init(2000*tpy, 2100*tpy);
		sim.simulate();

		vector<double> ba = sim.props.species.basal_area_vec;
		for (auto& b : ba) b*=1e4;
		cout << setprecision(10) << "Basal areas [m2/Ha]: " << ba << '\n';

		// err = is_equal(ba, {1.568911254, 7.061621143, 27.19173438});  // this was the output at yearly step_to when using 1d libpspm, reproduced in 2d
		// err = is_equal(ba, {1.569239985, 7.069009343, 27.21039394});  // this is the output  at yearly step_to after ODE bugfix in libpspm @8d0eb68
		// err = is_equal(ba, {1.569239985, 7.069009343, 27.21039394});  // this is the output  at yearly step_to after ODE bugfix in libpspm @118d623
		// err = is_equal(ba, {1.569460052, 7.060768412, 27.17424297});  // this is the output  at BIWEEKLY step_to. This difference is due to insertion of an extra point in R0 moving averager at integer t in yearly stepping - this doesnt happen at biweekly stepping, which seems more correct. This is the best level of investigattion I can do for now, so I'm going to take this as the expected test result and proceed with other stuff.
		// err = is_equal(ba, {1.569063555, 7.059273716, 27.16867421});  // this is the output  at BIWEEKLY step_to after moving climate update to before step. Difference is because now climate is not updated in after_step() before computation of seed rain
		// err = is_equal(ba, {1.570225188, 7.046974989, 27.14456976});  // this is the output  at BIWEEKLY step_to after new climate stream impl. Difference is because the two versions skip different sets of climate indices due to floating point comparisons.
		// err = is_equal(ba, {1.579158228, 7.140990396, 26.88772544});  // this is the output  at BIWEEKLY step_to after upgrade to latest version of phydro @6fc30d6. Difference is due to new temperature dependencies in phydro
		// err = is_equal(ba, {1.723965758, 7.637388237, 26.86261917});  // this is the output  at BIWEEKLY step_to after imlepmenting instantaneous version of phydro @ea5b867. 
		// err = is_equal(ba, {1.724216996, 7.638630557, 26.87059782});  // this is the output  at BIWEEKLY step_to after setting acclimation timescale to 7 days.
		// err = is_equal(ba, {1.723887226, 7.637860375, 26.86129712});  // this is the output  at BIWEEKLY step_to after adding 1-e6 to t in climate reading, and setting acclim_tc as t_growth in inst model.
		// err = is_equal(ba, {1.723740236, 7.642284321, 26.8675819});   // this is the output  at BIWEEKLY step_to after moving calc_seed_rain_r0.out of afterStep() and after step_to(). Difference is because there are still some tiny steps. 
		// err = is_equal(ba, {1.723090486, 7.639859488, 26.85158121});  // this is the output  at BIWEEKLY step_to after using actual p88/p50 ratio, instead of 3.01, in coordinateTraits()
		// err = is_equal(ba, {1.651351438, 7.370938402, 26.55374903});  // this is the output  at BIWEEKLY step_to after using daylength of 0.5 in phydro
		// err = is_equal(ba, {1.651263056, 7.37018801, 26.55223556});  // this is the output  at BIWEEKLY step_to after bugfix in unit conversion in assimilation.tpp
		err = is_equal(ba, {1.651105762, 7.36321835, 26.55250009});  // this is the output  at BIWEEKLY step_to after adding this line in simulate_t(): if (fabs(t - S.current_time) < 1e-6) continue;

		sim.close();
	}

	return err;

	// for (int i=0; i<1; ++i){
	// 	Patch sim(pfile);
	// 	sim.expt_dir = sim.expt_dir + "_614ppm";
	// 	sim.E.clim_inst.co2 = 614;
	// 	sim.init(1000, 5000);
	// 	sim.simulate();
	// 	sim.close();
	// }


	// // 1. eCO2 run for 2 different zetas
	// vector<double> zeta = {0.2, 0.08};

	// for (int i=0; i<1; ++i){
	// 	Patch sim(pfile);
	// 	sim.traits0.zeta = zeta[i];
	// 	sim.expt_dir = "HIST_ELE_zeta_" + to_string(zeta[i]);
	// 	sim.init(-1000, 5000);
	// 	sim.simulate();
	// 	sim.close();
	// }


	// // Test effect of zeta
	// vector<double> zeta = myseq(0.03, 0.3, 13);

	// std::for_each(
	// 	std::execution::par_unseq,
	// 	zeta.begin(),
	// 	zeta.end(),
	// 	[](double zz){
	// 		Patch sim("tests/params/p.ini");
	// 		sim.traits0.zeta = zz;
	// 		sim.expt_dir = "HIST_zeta_" + to_string(zz);
	// 		sim.init(1000, 1500);
	// 		sim.simulate();
	// 		sim.close();
	// 	});




	// // Effect of CO2
	// vector<double> co2_vec = myseq(360, 600, 12);
	// for (auto cc : co2_vec){
	// 	Patch sim("tests/params/p.ini");
	// 	sim.expt_dir = "scan_co2_" + to_string(cc);
	// 	sim.E.clim_inst.co2 = cc;
	// 	sim.init(1000, 1500);
	// 	sim.simulate();
	// 	sim.close();
	// }


	// Ensure that 1 mo timestep is good enough... 
	// {
	// 	Patch sim("tests/params/p.ini");
	// 	sim.timestep = 1.0/12.0/2.0;
	// 	sim.expt_dir = "BASE_dt_2wk";
	// 	sim.init(1000, 2000);
	// 	sim.simulate();
	// 	sim.close();
	// }

	// {
	// 	Patch sim("tests/params/p.ini");
	// 	sim.timestep = 1.0/12.0;
	// 	sim.expt_dir = "BASE_dt_4wk";
	// 	sim.init(1000, 2000);
	// 	sim.simulate();
	// 	sim.close();
	// }



	// Test effect of eCO2 under different zeta
	// {
	// 	Patch sim("tests/params/p_ele_base.ini");
	// 	sim.init(1000, 2000);
	// 	sim.simulate();
	// 	sim.close();
	// }

	// {
	// 	Patch sim("tests/params/p_ele_hi.ini");
	// 	sim.init(1000, 3000);
	// 	sim.simulate();
	// 	sim.close();
	// }

	// {
	// 	Patch sim("tests/params/p_ele_low.ini");
	// 	sim.init(1000, 3000);
	// 	sim.simulate();
	// 	sim.close();
	// }


}


