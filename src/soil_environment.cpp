#include "soil_environment.h"
#include <io_utils.h>

#include <numeric>
namespace pfate{
namespace env{

inline void print_vv(vector<vector<double>>& result){
	cout << setprecision(6);
	for (int d=0; d<6; ++d){
		for (auto& var : result) cout << var[d] << "\t";
		cout << '\n';
	}
	for (int d=359; d<result[0].size(); ++d){
		for (auto& var : result) cout << var[d] << "\t";
		cout << '\n';
	}	
}

SoilEnvironment::SoilEnvironment() : splash(0,0), solar(0,0) {
}

void SoilEnvironment::init(){
	splash = SPLASH(par_spl.lat, par_spl.elev);
	solar = SOLAR(par_spl.lat, par_spl.elev);
}

void SoilEnvironment::spinup(int y, std::vector<double>& sw_in, std::vector<double>& tair, std::vector<double>& pn, std::vector<double>& snow){
	
	double precip_total = std::accumulate(pn.begin(), pn.end(), 0.0);
	double snow_total   = std::accumulate(snow.begin(), snow.end(), 0.0);
	double Pinit = precip_total + snow_total;
	// 1. run spinup with AI=0 (garbage), to get an initial estimate of AI
	// result: wn_vec (0), snow_vec (1), qin_prev_vec (2), tdrain_vec (3), ro_vec (4), nds_prev_vec (5), pet_vec (6)
	vector<vector<double>> result = splash.spin_up_cpp(sw_in.size(), y, sw_in, tair, pn, 
                                                       par_spl.slope, par_spl.asp, snow, par_spl.soil_info);

	double pet_total = std::accumulate(result[6].begin(), result[6].end(), 0.0); 

	double AI_init = pet_total/Pinit;

	cout << "Initial spinup output:\n";
	print_vv(result);
	cout << "Initial AI = " << AI_init << '\n';

	par_spl.soil_info[12] = AI_init;
	result = splash.spin_up_cpp(sw_in.size(), y, sw_in, tair, pn, 
								par_spl.slope, par_spl.asp, snow, par_spl.soil_info);

	cout << "Final spinup output:\n";
	print_vv(result);

	// Save last state from spinup so that we can continue for main run
	state.wn = result[0].back();
	state.swe = result[1].back();
	state.qin = result[2].back();
	state.td = result[3].back();
	state.nd = result[5].back();
}

void SoilEnvironment::update_radiation(int doy, int y, double sw_in, double tair, double pn, double snow){
	solar.calculate_daily_fluxes(doy, y, sw_in, tair, par_spl.slope, par_spl.asp, snow, state.nd+1, dsoil.stress_factor);
	//                                                                                           ^ FIXME: This is just a hacky fix to temporarily ensure that nd seen by solar calc in SPLASH is the same as that seen by it at start of step. But this should be fixed INSIDE SPLASH, see Q for david              
}

void SoilEnvironment::water_balance_splash(int doy, int y, double sw_in, double tair, double pn, double snow, double plant_uptake){
	dvap = splash.run_one_day(doy, y, state.wn, sw_in, tair, pn, dsoil, par_spl.slope, par_spl.asp, state.swe, snow, par_spl.soil_info, state.qin, state.td, state.nd, plant_uptake);

	// cout << "td = " << dsoil.tdr << "\n";
	// FIXME: A temporary hack, needs to be fixed in SPLASH
	if (isnan(dsoil.tdr)) dsoil.tdr = 0;

	state.wn = dsoil.sm;
	state.swe = dsoil.swe;
	state.qin = dsoil.sqout;
	state.td = dsoil.tdr;
	state.nd = dsoil.nd;
}

void SoilEnvironment::water_balance(double dt, double precip, double pe_soil){
	double runoff = 0;

	double f = swc/par.rzwsc;

	double ae_soil = pe_soil * f;

	swc += (precip - ae_soil)*dt;

	if (swc < 0) swc = 0;

	if (swc > par.rzwsc){
		runoff = swc - par.rzwsc;
		swc = par.rzwsc;
	}
}

double SoilEnvironment::get_swp(){
	return pow(par.rzwsc/par.sstar, -par.b) - pow(swc/par.sstar, -par.b);
}

// Change log:
// v2: add dsoil
void SoilEnvironment::save(std::ostream& fout){
	fout << "SoilEnvironment::v2" << '\n';

	// output SoilPar
	fout << std::make_tuple(
		  par.b
		, par.rzwsc
		, par.sstar)
		<< '\n';

	// output SplashPars
	fout << std::make_tuple(
		    par_spl.slope
		  , par_spl.asp
		  , par_spl.lat
		  , par_spl.elev)
		<< '\n';
	
	fout << par_spl.soil_info << '\n';

	// output state
	fout << std::make_tuple(
		  state.wn
		, state.swe
		, state.qin
		, state.td
		, state.nd)
		<< '\n';

	// output dsoil (some of these are already in state, but thats ok)
	fout << std::make_tuple(
		  dsoil.sm
		, dsoil.ro
		, dsoil.swe
		, dsoil.bflow
		, dsoil.sqout
		, dsoil.tdr
		, dsoil.nd
		, dsoil.pet
		, dsoil.stress_factor
		, dsoil.psi_m
		, dsoil.sm_vv)
		<< '\n';


	fout << std::make_tuple(
		swc)
		<< '\n';

}


void SoilEnvironment::restore(std::istream& fin){
	std::string s; fin >> s; // discard version number
	assert(s == "SoilEnvironment::v2");

	fin >> par.b
		>> par.rzwsc
		>> par.sstar;

	fin >> par_spl.slope
		>> par_spl.asp
		>> par_spl.lat
		>> par_spl.elev;

	fin >> par_spl.soil_info;

	fin >> state.wn
		>> state.swe
		>> state.qin
		>> state.td
		>> state.nd;

	fin >> dsoil.sm
		>> dsoil.ro
		>> dsoil.swe
		>> dsoil.bflow
		>> dsoil.sqout
		>> dsoil.tdr
		>> dsoil.nd
		>> dsoil.pet
		>> dsoil.stress_factor
		>> dsoil.psi_m
		>> dsoil.sm_vv;

	fin >> swc;

}


} // env
} // namespace pfate
