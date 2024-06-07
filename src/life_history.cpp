#include "life_history.h"
#include <io_utils.h>
using namespace std;

namespace pfate{

ErgodicEnvironment::ErgodicEnvironment() : LightEnvironment(), Climate(){
	z_star = {15, 10, 5, 0};
	canopy_openness = {1, exp(-0.5 * 1.8), exp(-0.5 * 3.5), exp(-0.5 * 5.5)};
}

void ErgodicEnvironment::print(double t){
	Climate::print_line(t);
	// LightEnvironment::print();
	cout << "z_star = " << z_star;
	cout << "canopy_openness = " << canopy_openness;
}

void ErgodicEnvironment::computeEnv(double t, Solver* sol, std::vector<double>::iterator S, std::vector<double>::iterator dSdt){
	// do nothing
}


LifeHistoryOptimizer::LifeHistoryOptimizer(std::string params_file){
	//paramsFile = params_file; // = "tests/params/p.ini";
	I.parse(params_file);

	ts.set_units(I.get_verbatim("time_unit"));

	traits0.init(I);
	par0.init(I);

	par0.set_tscale(ts.get_tscale()); // default time unit is year

	c_stream.i_metFile = ""; //"tests/data/MetData_AmzFACE_Monthly_2000_2015_PlantFATE.csv";
	c_stream.a_metFile = ""; //"tests/data/MetData_AmzFACE_Monthly_2000_2015_PlantFATE.csv";
	c_stream.co2File = ""; //"tests/data/CO2_AMB_AmzFACE2000_2100.csv";

}

void LifeHistoryOptimizer::set_i_metFile(std::string file){
	c_stream.i_metFile = file;
	c_stream.update_i_met = (file == "") ? false : true;
}

void LifeHistoryOptimizer::set_a_metFile(std::string file){
	c_stream.a_metFile = file;
	c_stream.update_a_met = (file == "") ? false : true;
}

void LifeHistoryOptimizer::set_co2File(std::string co2file){
	c_stream.co2File = co2file;
	c_stream.update_co2 = (co2file == "") ? false : true;
}

void LifeHistoryOptimizer::init_co2(double _co2){
	C.init_co2(_co2);
}

void LifeHistoryOptimizer::init(){
	rep = 0;
	litter_pool = 0;
	seeds = 0;
	prod = 0;

	C.n_layers = C.z_star.size() - 1;

	C.set_elevation(0);
	C.set_acclim_timescale(7);
	c_stream.init();

	// We are tracking the life-cycle of a seed: how many seeds does a single seed produce (having gone through dispersal, germination, and plant life stages)
	P = plant::Plant();
	// P.initFromFile(params_file);
	P.init(par0, traits0);

	P.geometry.set_lai(P.par.lai0);
	P.set_size(0.01);
	// Simulation below starts at seedling stage. So account for survival until seedling stage
	P.state.mortality = -log(P.p_survival_dispersal(C) * P.p_survival_germination(C)); // p{fresh seed is still alive after germination} = p{it survives dispersal}*p{it survives germination}

	// double total_prod = P.get_biomass();
	// cout << "Starting biomass = " << total_prod << "\n";
	// cout << "Mortality until seedling stage = " << P.state.mortality << "\n";

}

vector<std::string> LifeHistoryOptimizer::get_header(){
	return {
		  "i"
		, "ppfd"
		, "assim_net"
		, "assim_gross"
		, "rl"
		, "rr"
		, "rs"
		, "tl"
		, "tr"
		, "dpsi"
		, "vcmax"
		, "transpiration"
		, "height"
		, "diameter"
		, "crown_area"
		, "lai"
		, "sapwood_fraction"
		, "leaf_mass"
		, "root_mass"
		, "stem_mass"
		, "coarse_root_mass"
		, "total_mass"
		, "total_rep"
		// , "seed_pool"
		// , "germinated"
		, "fitness"
		, "total_prod"
		, "litter_mass"
		, "mortality"
		, "mortality_inst"
		, "mortrate_0"
		, "mortrate_growth"
		, "mortrate_d"
		, "mortrate_hyd"
		, "leaf_lifespan"
		, "fineroot_lifespan"
	};
}

void LifeHistoryOptimizer::printHeader(ostream& lfout){
	vector<std::string> s = get_header();
	for (auto& vv : s) lfout << vv << "\t";
	lfout << '\n';
}

vector<double> LifeHistoryOptimizer::get_state(double t){
	return {
		  ts.to_julian(t)
		, C.clim_inst.ppfd
		, P.assimilator.plant_assim.npp
		, P.assimilator.plant_assim.gpp
		, P.assimilator.plant_assim.rleaf
		, P.assimilator.plant_assim.rroot
		, P.assimilator.plant_assim.rstem
		, P.assimilator.plant_assim.tleaf
		, P.assimilator.plant_assim.troot
		, P.assimilator.plant_assim.dpsi_avg
		, P.assimilator.plant_assim.vcmax_avg
		, P.assimilator.plant_assim.trans
		, P.geometry.height
		, P.geometry.diameter
		, P.geometry.crown_area
		, P.geometry.lai
		, P.geometry.sapwood_fraction
		, P.geometry.leaf_mass(P.traits)
		, P.geometry.root_mass(P.traits)
		, P.geometry.stem_mass(P.traits)
		, P.geometry.coarse_root_mass(P.traits)
		, P.get_biomass()
		, rep
	//  , P.state.seed_pool
	//  , germinated
		, seeds
		, prod
		, litter_pool
		, P.state.mortality
		, P.rates.dmort_dt
		, P.mort.mu_0
		, P.mort.mu_growth
		, P.mort.mu_d
		, P.mort.mu_hyd
		, 1 / P.assimilator.kappa_l
		, 1 / P.assimilator.kappa_r
	};
}

void LifeHistoryOptimizer::printState(double t, ostream& lfout){
	vector<double> v = get_state(t);
	for (auto& vv : v) lfout << vv << "\t";
	lfout << '\n';
}

// void LifeHistoryOptimizer::set_traits(std::vector<double> tvec){
// 	P.set_evolvableTraits(tvec);
// }


// std::vector<double> LifeHistoryOptimizer::get_traits(){
// 	return P.get_evolvableTraits();
// }


void LifeHistoryOptimizer::printMeta(){
	P.print();
	C.print(0);
}


void LifeHistoryOptimizer::set_state(vector<double>::iterator it){
	P.geometry.set_lai(*it++);
	P.set_size(*it++);
	prod = *it++;
	litter_pool = *it++;
	rep = *it++;
//		P.state.seed_pool = *it++;
	seeds = *it++;
	P.state.mortality = *it++;
}

void LifeHistoryOptimizer::get_rates(vector<double>::iterator it){
	*it++ = P.rates.dlai_dt;       // lai growth rate
	*it++ = P.rates.dsize_dt;    // size (diameter) growth rate
	*it++ = P.bp.dmass_dt_tot;	   // biomass production rate
	*it++ = P.bp.dmass_dt_lit;  // litter biomass growth rate
	*it++ = P.bp.dmass_dt_rep; //(1-fg)dBdt;  // reproduction biomass growth rate
//		*it++ = P.rates.dseeds_dt_pool;
	*it++ = P.rates.dseeds_dt;
	*it++ = P.rates.dmort_dt;
}


void LifeHistoryOptimizer::update_climate(double julian_time){
	c_stream.updateClimate(julian_time, C);
}

void LifeHistoryOptimizer::grow_for_dt(double t, double dt){

	auto derivs = [this](double t, std::vector<double>& S, std::vector<double>& dSdt){
		//if (fabs(t - 2050) < 1e-5) 
		double j_clim = ts.to_julian(t);
		update_climate(j_clim);
		C.set_forcing_acclim(j_clim, C.clim_midday);

		// C.Climate::print(t);
		set_state(S.begin());
		P.calc_demographic_rates(C, t);

		// Override Plant-FATE fecundity calculations 
		// We need to explicitly include plant mortality here for fitness calcs
		double fec = P.fecundity_rate(P.bp.dmass_dt_rep, C);
		P.rates.dseeds_dt =  fec * exp(-P.state.mortality);  // Fresh seeds produced = fecundity rate * p{plant is alive}
		// P.rates.dseeds_dt_germ =   P.state.seed_pool/P.par.ll_seed;   // seeds that leave seed pool proceed for germincation

		get_rates(dSdt.begin());
	};

	std::vector<double> S = {P.geometry.lai, P.geometry.get_size(), prod, litter_pool, rep, seeds, P.state.mortality};
	RK4(t, dt, S, derivs);
	//Euler(t, dt, S, derivs);
	set_state(S.begin());
}


double LifeHistoryOptimizer::calcFitness(){
	// lho_set_traits(tvec);
	for (double t=2000; t <= 2500; t=t + dt){
		grow_for_dt(t, dt);
	}
	return seeds;
}

} // namespace pfate
