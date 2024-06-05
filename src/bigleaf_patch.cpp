#include "bigleaf_patch.h"
#include <io_utils.h>
#include <SOLAR.h>
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

	// init soil params and site info for splash
	soil_env.par_spl.soil_info = I.get_vector<double>("soil_info");
	soil_env.par_spl.slope = I.get<double>("slope");
	soil_env.par_spl.asp = 180-I.get<double>("aspect");
	soil_env.par_spl.lat = I.get<double>("latitude");
	soil_env.par_spl.elev = I.get<double>("elevation");
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

void BigLeafPatch::init(double _t0, double _tf){
	t0 = _t0; 
	tf = _tf;
	forcing.set_elevation(0);
	forcing.set_acclim_timescale(7);
	climate_stream.init();
	soil_env.swc = 50;
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
		, "precip_acclim"
		, "tc_inst"
		, "ppfd_inst"
		, "rn_inst"
		, "vpd_inst"
		, "co2_inst"
		, "elv_inst"
		, "swp_inst"
		, "precip_inst"
		, "swc"
		, "assim_gross"
		, "dpsi"
		, "vcmax"
		, "jmax"
		, "gs"
		, "transpiration"
		, "latent_energy"
		, "le_wet_soil"
		, "soil_mois"
		, "swp"
		, "aet_splash"
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
		, forcing.clim_acclim.precip
		, forcing.clim_inst.tc
		, forcing.clim_inst.ppfd
		, forcing.clim_inst.rn
		, forcing.clim_inst.vpd
		, forcing.clim_inst.co2
		, forcing.clim_inst.elv
		, forcing.clim_inst.swp
		, forcing.clim_inst.precip
		, soil_env.swc
		, phydro_out.a
		, phydro_out.dpsi
		, phydro_out.vcmax
		, phydro_out.jmax
		, phydro_out.gs
		, phydro_out.e
		, phydro_out.le
		, soil_env.dvap.aet*(phydro_out.le/phydro_out.e) // phydro_out.le_s_wet
		, soil_env.state.wn
		, soil_env.dsoil.psi_m
		, soil_env.dvap.aet/(1-fapar)
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

void BigLeafPatch::simulate(){

	ofstream fout(outfile.c_str());
	printHeader(fout);
	cout << "Main run: " << t0 << " --> " << tf << " (dt = " << timestep << ")\n";
	for (double t=t0+timestep*0.5; t <= tf + 1e-6; t=t + timestep) {  // 1e-6 ensures that last timestep to reach yf is actually executed
		auto time_pt = flare::julian_to_date(ts.to_julian(t));
		int year = time_pt.tm_year+1900;
		int doy = time_pt.tm_yday+1;
		cout << flare::julian_to_datestring(ts.to_julian(t)) << " (" << year << "." << doy << "): "; 

		// read forcing inputs
		// std::cout << "update Env (explicit)... t = " << S.current_time << ":\n";
		update_climate(ts.to_julian(t) + 1e-6); // The 1e-6 is to ensure that when t coincides exactly with time in climate file, we ensure that the value in climate file is read by asking for a slightly higher t
		// forcing.print_line(t);

		// start of step stress factor
		cout << "start of step stress factor = " << soil_env.dsoil.stress_factor << '\n';

		// calc radiation components
		// FIXME: update nd here
		soil_env.update_radiation(doy, year, forcing.clim_inst.ppfd/2.04, forcing.clim_inst.tc, forcing.clim_inst.precip, 0);
		srad sol = soil_env.solar.get_vals();
	    cout << "Rnl in before func: " << sol.rnl << "W/m^2\n";
	    cout << "ppfd_in/ppfd_net = " << forcing.clim_inst.ppfd << " / " << sol.ppfd_d*1e6/86400 << '\n';

		forcing.set_forcing_acclim(ts.to_julian(t) + 1e-6, forcing.clim_midday);

		// forcing.clim_inst.swp = soil_env.get_swp();
		// forcing.clim_acclim.swp = soil_env.get_swp();
		forcing.clim_inst.swp = soil_env.dsoil.psi_m;
		forcing.clim_acclim.swp = soil_env.dsoil.psi_m;

		// photosynthesis
		phydro_out = bigleaf_assimilator.leaf_assimilation_rate(1, fapar, forcing, par0, traits0);
		
		// Convert units 
		phydro_out.a     *= (86400 * 1e-6 * 12);  // umol co2/m2/s ----> umol co2/m2/day --> mol co2/m2/day --->  gC /m2/day
		phydro_out.e     *= (86400 * 0.018015 * 1e-3 * 1000);   // mol h2o/m2/s ---> mol h2o/m2/day ---> kg h2o/m2/day ---> m3 /m2/day ---> mm/day   
		phydro_out.le    *= 86400;  // W m-2 = J m-2 s-1 ---> J m-2 day-1
		phydro_out.le_s_wet  *= 86400;  // W m-2 = J m-2 s-1 ---> J m-2 day-1

		// water balance
		double m = par0.days_per_tunit; // multiplier to convert from day-1 to t_unit-1
		soil_env.water_balance(timestep, forcing.clim_inst.precip*m, phydro_out.e*m);

		soil_env.water_balance_splash(doy, year, forcing.clim_inst.ppfd/2.04*(1-fapar), forcing.clim_inst.tc, forcing.clim_inst.precip, 0, phydro_out.e);
		
		cout << soil_env.state.wn << " " << soil_env.dsoil.ro << " " << soil_env.dsoil.ro << " " << soil_env.state.nd << " " << soil_env.dsoil.stress_factor << " " << soil_env.dsoil.psi_m << "\n";

		printState(t, fout);
	}
	fout.close();
}


void BigLeafPatch::spinup(){

	ofstream fout(outfile.c_str());
	printHeader(fout);
	vector<double> sw_in, tair, precip, snowfall;

	// 1. Collect daily 24-hr mean forcing variables for 1 year for soil moisture spinup
	double dt_spin = 1/par0.days_per_tunit; // 1 day timstep for spinup
	double tend = t0 + 1/par0.years_per_tunit_avg; // 1 year of data for spinup run
	cout << "Spinup: " << t0 << " --> " << tend << " (dt = " << dt_spin << ")\n";
	auto tp = flare::julian_to_date(ts.to_julian(t0+dt_spin));
	int y = tp.tm_year+1900;
	int m = tp.tm_mon+1;
	int d = tp.tm_mday;
	if (m > 1 && d > 1) throw std::runtime_error("Spin up should start on the first day of the year");

	cout << "Spinup: start date = " << y << "-" << m << "-" << d << "\n";
	for (double t=t0+dt_spin*0.5; t <= tend + 1e-6; t=t + dt_spin) {  // 1e-6 ensures that last timestep to reach yf is actually executed
		// read forcing inputs
		update_climate(ts.to_julian(t) + 1e-6); // The 1e-6 is to ensure that when t coincides exactly with time in climate file, we ensure that the value in climate file is read by asking for a slightly higher t
		// forcing.print_line(t);

		sw_in.push_back(forcing.clim_inst.ppfd/2.04);
		tair.push_back(forcing.clim_inst.tc);
		precip.push_back(forcing.clim_inst.precip);
		snowfall.push_back(0);

		cout << flare::julian_to_datestring(ts.to_julian(t)) << "\t" 
		     << sw_in.back() << "\t" 
		     << tair.back() << "\t" 
		     << precip.back() << "\t" 
		     << snowfall.back() << "\n"; 
	}

	// 2. Run spinup
	soil_env.spinup(y, sw_in, tair, precip, snowfall);

	// 3. Since main run will start from back in time, we need to reset the moving averager for acclim forcing
	forcing.reset_forcing_acclim();
}


} // namespace pfate
