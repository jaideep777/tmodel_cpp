#include "soil_environment.h"
#include <SPLASH.h>

#include <numeric>
namespace pfate{
namespace env{

inline void print_vv(vector<vector<double>>& result){
	for (int d=0; d<6; ++d){
		for (auto& var : result) cout << var[d] << "\t";
		cout << '\n';
	}
	for (int d=359; d<365; ++d){
		for (auto& var : result) cout << var[d] << "\t";
		cout << '\n';
	}	
}

void SoilEnvironment::spinup(int y, std::vector<double>& sw_in, std::vector<double>& tair, std::vector<double>& pn, std::vector<double>& snow){
	SPLASH splash(par_spl.lat, par_spl.elev);
	
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



} // env
} // namespace pfate
