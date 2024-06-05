#include "climate_stream.h"

namespace pfate{
namespace env{

void ClimateStream::init(){
	if (update_co2){
		co2_stream.set_tname("Year");
		co2_stream.periodic = false;
		co2_stream.centered_t = false;
		co2_stream.open({co2File}, "years CE");
	}
	if (update_i_met){
		i_met_stream.set_tname("Decimal_year");
		i_met_stream.periodic = true;
		i_met_stream.centered_t = false;
		i_met_stream.open({i_metFile}, "years CE");
	}
	if (update_a_met){
		a_met_stream.set_tname("Decimal_year");
		a_met_stream.periodic = true;
		a_met_stream.centered_t = false;
		a_met_stream.open({a_metFile}, "years CE");
	}
}


void ClimateStream::updateClimate(double julian_time, Climate& C){
	if (update_co2){
		co2_stream.advance_to_time(julian_time);
		// std::cout << co2_stream.current_row << std::endl;
		C.clim_inst.co2   = as<double>(co2_stream.current_row[1]); // co2 is in index 1
		C.clim_acclim.co2 = as<double>(co2_stream.current_row[1]); // co2 is in index 1
	}
	if (update_i_met){
		i_met_stream.advance_to_time(julian_time);
		// std::cout << i_met_stream.current_row << std::endl;
		C.clim_inst.tc   = as<double>(i_met_stream.current_row[3]);
		C.clim_inst.vpd  = as<double>(i_met_stream.current_row[4]) * 100; // convert hPa to Pa
		C.clim_inst.ppfd = as<double>(i_met_stream.current_row[5]);      // ppfd
		C.clim_inst.swp  = as<double>(i_met_stream.current_row[7]) * (-1); // convert -MPa to MPa
		C.clim_inst.rn   = C.clim_inst.ppfd/2; // set to half of ppfd for now, eventually use solar model to calculate
		if (use_precip_data) C.clim_inst.precip = as<double>(i_met_stream.current_row[8]); // precipitation (optional)
	}
	if (update_a_met){
		a_met_stream.advance_to_time(julian_time);
		// std::cout << a_met_stream.current_row << std::endl;
		Clim cnew = C.clim_acclim;
		cnew.tc   = as<double>(a_met_stream.current_row[3]);
		cnew.vpd  = as<double>(a_met_stream.current_row[4]) * 100; // convert hPa to Pa
		cnew.ppfd = as<double>(a_met_stream.current_row[6]);      // max ppfd
		cnew.swp  = as<double>(a_met_stream.current_row[7]) * (-1); // convert -MPa to MPa
		cnew.rn   = cnew.ppfd/2; // set to half of ppfd for now, eventually use solar model to calculate
		cnew.precip = 0; // precipitation is never used in acclim data
		C.clim_midday = cnew;
		// C.set_forcing_acclim(julian_time, cnew);
	}
	else {
		// TODO. Replace with a proper function that computes acclim forcing from inst
		Clim cnew = C.clim_acclim;
		cnew = C.clim_inst;
		cnew.ppfd = C.clim_inst.ppfd * 4;
		cnew.rn   = cnew.ppfd/2; // set to half of ppfd for now, eventually use solar model to calculate
		C.clim_midday = cnew;
		// C.set_forcing_acclim(julian_time, cnew);
	}

}

} // namespace env
} // namespace pfate
