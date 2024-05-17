#include "bigleaf_patch.h"
#include <io_utils.h>
using namespace std;

namespace pfate{


BigLeafPatch::BigLeafPatch(std::string params_file){
	//paramsFile = params_file; // = "tests/params/p.ini";
	I.parse(params_file);

	traits0.init(I);
	par0.init(I);

	timestep = I.get<double>("timestep");  // timestep	
	ts.set_units(I.get_verbatim("time_unit"));
	par0.set_tscale(ts.get_tscale()); // default time unit is year

	climate_stream.i_metFile = I.get<string>("i_metFile");
	climate_stream.a_metFile = I.get<string>("a_metFile");
	climate_stream.co2File   = I.get<string>("co2File");
	climate_stream.update_i_met = (climate_stream.i_metFile == "null") ? false : true;
	climate_stream.update_a_met = (climate_stream.a_metFile == "null") ? false : true;
	climate_stream.update_co2   = (climate_stream.co2File == "null") ? false : true;

	outfile = I.get<string>("outDir") + "/" 
	        + I.get<string>("exptName") + "/" 
			+ "bigleaf_" + I.get<std::string>("emgProps");

	traits0.init(I);
	traits0.coordinate(); 

	par0.init(I);

	traits0.print();
	par0.print();
}

void BigLeafPatch::set_i_metFile(std::string file){
	climate_stream.i_metFile = file;
	climate_stream.update_i_met = (file == "") ? false : true;
}

void BigLeafPatch::set_a_metFile(std::string file){
	climate_stream.a_metFile = file;
	climate_stream.update_a_met = (file == "") ? false : true;
}

void BigLeafPatch::set_co2File(std::string co2file){
	climate_stream.co2File = co2file;
	climate_stream.update_co2 = (co2file == "") ? false : true;
}

void BigLeafPatch::init_co2(double _co2){
	forcing.init_co2(_co2);
}

void BigLeafPatch::set_fapar(double _fapar){
	fapar = _fapar;
}

void BigLeafPatch::init(){

	forcing.set_elevation(0);
	forcing.set_acclim_timescale(7);
	climate_stream.init();

}

vector<std::string> BigLeafPatch::get_header(){
	return {
		  "i"
		, "tc_acclim"
		, "ppfd_acclim"
		, "rn_acclim"
		, "vpd_acclim"
		, "co2_acclim"
		, "elv_acclim"
		, "swp_acclim"
		, "tc_inst"
		, "ppfd_inst"
		, "rn_inst"
		, "vpd_inst"
		, "co2_inst"
		, "elv_inst"
		, "swp_inst"
		, "assim_gross"
		, "dpsi"
		, "vcmax"
		, "jmax"
		, "gs"
		, "transpiration"
	};
}

void BigLeafPatch::printHeader(ostream& lfout){
	vector<std::string> s = get_header();
	for (auto& vv : s) lfout << vv << "\t";
	lfout << '\n';
}

vector<double> BigLeafPatch::get_state(double t){
	return {
		  forcing.clim_acclim.tc
		, forcing.clim_acclim.ppfd
		, forcing.clim_acclim.rn
		, forcing.clim_acclim.vpd
		, forcing.clim_acclim.co2
		, forcing.clim_acclim.elv
		, forcing.clim_acclim.swp
		, forcing.clim_inst.tc
		, forcing.clim_inst.ppfd
		, forcing.clim_inst.rn
		, forcing.clim_inst.vpd
		, forcing.clim_inst.co2
		, forcing.clim_inst.elv
		, forcing.clim_inst.swp
		, phydro_out.a
		, phydro_out.dpsi
		, phydro_out.vcmax
		, phydro_out.jmax
		, phydro_out.gs
		, phydro_out.e
	};
}

void BigLeafPatch::printState(double t, ostream& lfout){
	vector<double> v = get_state(t);
	lfout << flare::date_to_string(flare::julian_to_date(ts.to_julian(t)), "%Y-%m-%d %H:%M:%S") << '\t';
	for (auto& vv : v) lfout << vv << "\t";
	lfout << '\n';
}

void BigLeafPatch::printMeta(){
	forcing.print(0);
}


void BigLeafPatch::update_climate(double julian_time){
	climate_stream.updateClimate(julian_time, forcing);
}

void BigLeafPatch::simulate(double t0, double tf){

	ofstream fout(outfile.c_str());
	printHeader(fout);
	for (double t=t0; t <= tf + 1e-6; t=t + timestep) {  // 1e-6 ensures that last timestep to reach yf is actually executed
		// read forcing inputs
		// std::cout << "update Env (explicit)... t = " << S.current_time << ":\n";
		update_climate(ts.to_julian(t) + 1e-6); // The 1e-6 is to ensure that when t coincides exactly with time in climate file, we ensure that the value in climate file is read by asking for a slightly higher t

		// photosynthesis
		phydro_out = bigleaf_assimilator.leaf_assimilation_rate(1, fapar, forcing, par0, traits0);
		
		// Convert units 
		phydro_out.a     *= (86400 * 1e-6 * 12);  // umol co2/m2/s ----> umol co2/m2/day --> mol co2/m2/day --->  gC /m2/day
		phydro_out.e     *= (86400 * 0.018015 * 1e-3 * 1000);   // mol h2o/m2/s ---> mol h2o/m2/day ---> kg h2o/m2/day ---> m3 /m2/day ---> mm/day   
		phydro_out.le    *= 86400;  // W m-2 = J m-2 s-1 ---> J m-2 day-1

		// water balance
		

		printState(t, fout);
	}
	fout.close();
}


} // namespace pfate
