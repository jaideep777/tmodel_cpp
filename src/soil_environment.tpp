namespace pfate{
namespace env{

template <class ForcingInput>
void SoilEnvironment::update_forcings(int doy, int year, ForcingInput& forcing){
	// TODO. calc snowfall from p_snow (Currently set to 0), and update precip to liquid rain
	double p_snow = 0; // probability of snowfall // TODO: replace with proper function
	double f_rain = 1; // fraction of precip falling as rain (TODO: replace with proper function)

	forcing.clim_inst.snow = forcing.clim_inst.precip*(1-f_rain);
	forcing.clim_inst.precip *= f_rain;

	// update snow age, as it's required by radiation calcs
	//    only a copy is updated here, the actual state variable will be updated during water balance calcs
	double nd_today = (forcing.clim_inst.snow > 0.0)? 0 : state.nd+1;

	// calc radiation components
	solar.calculate_daily_fluxes(doy, year, forcing.clim_inst.ppfd/2.04, forcing.clim_inst.tc, par_spl.slope, par_spl.asp, forcing.clim_inst.snow, nd_today, dsoil.stress_factor);
	srad sol = solar.get_vals();
	// cout << "Rnl in before func: " << sol.rnl << "W/m^2\n";
	// cout << "ppfd_in/ppfd_net = " << forcing.clim_inst.ppfd << " / " << sol.ppfd_d*1e6/86400 << '\n';

	// Set radiation and swp from SPLASH 
	double rn_midday_by_24hr = forcing.clim_midday.ppfd/forcing.clim_inst.ppfd;
	forcing.clim_inst.ppfd = sol.ppfd_d*1e6/86400; // convert mol m-2 day ---> umol m-2 s-1
	forcing.clim_midday.ppfd = forcing.clim_inst.ppfd*rn_midday_by_24hr;

	forcing.clim_inst.swp = dsoil.psi_m;
	forcing.clim_midday.swp = dsoil.psi_m;

	// cout << "swp = " << forcing.clim_inst.swp << " / " << forcing.clim_acclim.swp << '\n';

}


} // env
} // namespace pfate
