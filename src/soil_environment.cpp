#include "soil_environment.h"

namespace pfate{
namespace env{

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
