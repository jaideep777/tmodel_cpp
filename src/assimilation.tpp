namespace plant{

inline void print_phydro(const phydro::PHydroResult& res, std::string s){
	std::cout << "phydro: " << s << '\n';
	std::cout 
		<< "   a = " << res.a << '\n'
		<< "   e = " << res.e << '\n'
		<< "   vcmax = " << res.vcmax << '\n'
		<< "   vcmax25 = " << res.vcmax25 << '\n'
		<< "   rd = " << res.rd << '\n'
		<< "   gs = " << res.gs << '\n'
		<< "   dpsi = " << res.dpsi << '\n';
}

// **
// ** Gross and Net Assimilation 
// **
template<class _Climate>
phydro::PHydroResult Assimilator::leaf_assimilation_rate(double fipar, double fapar, _Climate& C, PlantParameters& par, PlantTraits& traits){

	// std::cout << "--- In LEAF ASIMILATION RATE ---" << std::endl;
	phydro::ParCost par_cost(par.alpha, par.gamma);
	phydro::ParPlant par_plant(traits.K_leaf, traits.p50_leaf, traits.b_leaf);
	phydro::ParControl par_control;

	par_control.gs_method = phydro::GS_APX;
	par_control.et_method = phydro::ET_DIFFUSION;

	double f_day_length = 0.5;

	double Iabs_acclim = fipar * C.clim_acclim.ppfd;
	double Iabs_day    = fipar * C.clim_inst.ppfd / f_day_length;
	double Iabs_24hr   = fipar * C.clim_inst.ppfd;
	// std::cout << "Irradiation inputs: acclim " << Iabs_acclim << " ; day " << Iabs_day << " ; 24 hr " << Iabs_24hr << std::endl;

	auto out_phydro_acclim = phydro::phydro_analytical(
		C.clim_acclim.tc,     // current temperature
		C.clim_acclim.tc,     // growth temperature
		Iabs_acclim,          // midday incident PAR [umol m-2 s-1]
		C.clim_acclim.rn,     // Net radiation [W m-2] (only used for LE calculations which we dont use) // FIXME. Should this be Rnl? See message to Beni
		C.clim_acclim.vpd,    // vpd [kPa]
		C.clim_acclim.co2,	  // co2 [ppm]
		C.clim_acclim.pa,     // surface pressure [Pa]
		fapar,                // fraction of absorbed PAR
		par.kphio,            // phi0 - quantum yield
		C.clim_acclim.swp,    // soil water potential [MPa]
		par.rd,               // ratio or dark respiration to vcmax
		C.clim_acclim.vwind,  // wind speed [m s-1], only used by PML, which we dont use, so set to global average of 3 m/s
		par_plant,            // plant hydraulic traits
		par_cost,             // cost params
		par_control           // configuration params for phydro
	);

	// print_phydro(out_phydro_acclim, "acclim"); //THIS WORKS

	// // ~~~~~
	// // Note: Inst calcs were being done in a hacky way before the instantaneous model was ready, as follows. 
	// // This was completely wrong!
	// auto photo_leaf1 = out_phydro_acclim;
	// // the factor f = 1.18 accounts for the non-linearity in the instantaneous sub-daily response in the P-hydro model
	// double f = 1.18*C.clim_inst.ppfd/C.clim_acclim.ppfd;
	// photo_leaf1.a *= f;
	// photo_leaf1.e *= f;
	// photo_leaf1.gs *= f;
	// // ~~~~~
	// print_phydro(photo_leaf1, "inst shortcut");

	// std::cout << "^^^ phydro instantaneous analytical input: ^^^" << std::endl;
	// std::cout << "out_phydro_acclim vcmax25\t" << out_phydro_acclim.vcmax25 << std::endl;
	// std::cout <<  "out_phydro_acclim jmax25\t" << out_phydro_acclim.jmax25 << std::endl;
	// std::cout << "climate instantaneous tc\t" << C.clim_inst.tc << std::endl;
	// std::cout << "climate acclimation tc\t" << C.clim_acclim.tc << std::endl;
	// std::cout << "Irradiation daily\t" << Iabs_day << std::endl;
	// std::cout << "climate instantaneous net radiation rn\t" << C.clim_inst.rn << std::endl;
	// std::cout << "climate instantaneous VPD tc\t" << C.clim_inst.vpd << std::endl;
	// std::cout << "climate instantaneous co2\t" << C.clim_inst.co2 << std::endl;
	// std::cout << "climate instantaneous pa\t" << C.clim_inst.pa << std::endl;
	// std::cout << "fraction absorbed PAR fapar\t" << fapar << std::endl;
	// std::cout << "quantum yield phi0\t" << par.kphio << std::endl;
	// std::cout << "climate instantaneous soil water potential SWP\t" << C.clim_inst.swp << std::endl;
	// std::cout << "ratio dark respiration to vcmax \t" << par.rd << std::endl;
	// std::cout << "climate instantaneous wind velocity vwind\t" << C.clim_inst.vwind << std::endl;


	auto photo_leaf = phydro::phydro_instantaneous_analytical(
		out_phydro_acclim.vcmax25, // acclimated vcmax25
		out_phydro_acclim.jmax25,  // acclimated jmax25
		C.clim_inst.tc,            // current temperature
		C.clim_acclim.tc,          // growth temperature
		Iabs_day,                  // daytime mean incident PAR [umol m-2 s-1]
		C.clim_inst.rn,            // mean net radiation [W m-2] (only used for LE calculations which we dont use)
		C.clim_inst.vpd,           // vpd [kPa]
		C.clim_inst.co2,	       // co2 [ppm]
		C.clim_inst.pa,            // surface pressure [Pa]
		fapar,                     // fraction of absorbed PAR
		par.kphio,                 // phi0 - quantum yield
		C.clim_inst.swp,           // soil water potential [MPa]
		par.rd,                    // ratio or dark respiration to vcmax
		C.clim_inst.vwind,         // wind speed [m s-1], only used by PML, which we dont use, so set to global average of 3 m/s
		par_plant,                 // plant hydraulic traits
		par_cost,                  // cost params
		par_control                // configuration params for phydro
	);

	photo_leaf.a        *= f_day_length;
	photo_leaf.e        *= f_day_length;
	photo_leaf.gs       *= f_day_length;
	// print_phydro(photo_leaf, "inst real 12 hr");


	// auto photo_leaf2 = phydro::phydro_instantaneous_analytical(
	// 	out_phydro_acclim.vcmax25, // acclimated vcmax25
	// 	out_phydro_acclim.jmax25,  // acclimated jmax25
	// 	C.clim_inst.tc,            // current temperature
	// 	C.clim_acclim.tc,          // growth temperature
	// 	Iabs_24hr,                 // 24-hr mean incident PAR [umol m-2 s-1]
	// 	C.clim_inst.rn,            // mean net radiation [W m-2] (only used for LE calculations which we dont use)
	// 	C.clim_inst.vpd,           // vpd [kPa]
	// 	C.clim_inst.co2,	       // co2 [ppm]
	// 	C.clim_inst.pa,            // surface pressure [Pa]
	// 	fapar,                     // fraction of absorbed PAR
	// 	par.kphio,                 // phi0 - quantum yield
	// 	C.clim_inst.swp,           // soil water potential [MPa]
	// 	par.rd,                    // ratio or dark respiration to vcmax
	// 	C.clim_inst.vwind,         // wind speed [m s-1], only used by PML, which we dont use, so set to global average of 3 m/s
	// 	par_plant,                 // plant hydraulic traits
	// 	par_cost,                  // cost params
	// 	par_control                // configuration params for phydro
	// );

	// print_phydro(photo_leaf2, "inst real 24 hr");

	// std::cout << "--- Finished LEAF ASIMILATION RATE ---" << std::endl;
	// std::cout << "--- --- --- --- ---" << std::endl;
	return photo_leaf;
}


template<class Env>
void  Assimilator::calc_plant_assimilation_rate(Env& env, PlantArchitecture* G, PlantParameters& par, PlantTraits& traits){
	//double GPP_plant = 0, Rl_plant = 0, dpsi_avg = 0;
	double fapar = 1 - exp(-par.k_light * G->lai);
	bool by_layer = false;

	plant_assim.gpp        = 0;
	plant_assim.rleaf      = 0;
	plant_assim.trans      = 0;
	plant_assim.dpsi_avg   = 0;
	plant_assim.vcmax_avg  = 0;
	plant_assim.vcmax25_avg = 0;
	plant_assim.mc_avg     = 0;
	plant_assim.gs_avg     = 0;
	plant_assim.c_open_avg = 0;

	double ca_cumm = 0;
	// std::cout << "--- PPA Assim begin ---" << "\n";
	for (int ilayer=0; ilayer <= env.n_layers; ++ilayer){ // for l in 1:layers{	
		double zst = env.z_star[ilayer];
		double ca_layer = G->crown_area_above(zst, traits) - ca_cumm;
		// std::cout << "h = " << G->height << ", z* = " << zst << ", I = " << env.canopy_openness[ilayer] << ", fapar = " << fapar << /*", A = " << (res.a + res.vcmax*par.rd) << " umol/m2/s x " <<*/ ", ca_layer = " << ca_layer << /*" m2 = " << (res.a + res.vcmax*par.rd) * ca_layer << ", vcmax = " << res.vcmax <<*/ "\n"; 
		// std::cout <<"by layer: " << by_layer << std::endl;
		// std::cout << "layer number: " << ilayer << std::endl;
		if (by_layer == true){
			auto res = leaf_assimilation_rate(env.canopy_openness[ilayer], fapar, env, par, traits);
			plant_assim.gpp        += (res.a + res.vcmax * par.rd) * ca_layer;
			plant_assim.rleaf      += (res.vcmax * par.rd) * ca_layer;
			plant_assim.trans      += res.e * ca_layer;
			plant_assim.dpsi_avg   += res.dpsi * ca_layer;
			plant_assim.vcmax_avg  += res.vcmax * ca_layer;
			plant_assim.gs_avg     += res.gs * ca_layer;
			plant_assim.vcmax25_avg += res.vcmax25 * ca_layer;
			plant_assim.mc_avg     += res.mc * ca_layer;
			// std::cout << "--- layer " << ilayer << "LEAF ASSIMILATION RATE begin ---"<< std::endl;
			// std::cout << "ca_layer " << ca_layer << std::endl;
			// std::cout << "a " << res.a << std::endl;
			// std::cout << "vcmax " << res.vcmax << std::endl;
			// std::cout << "rd " << res.rd << std::endl;
			// std::cout << "e " << res.e << std::endl;
			// std::cout << "dpsi " << res.dpsi << std::endl;
			// std::cout << "gs " << res.gs << std::endl;
			// std::cout << "vcmax25 " << res.vcmax25 << std::endl;
			// std::cout << "mc " << res.mc << std::endl;
		}

		plant_assim.c_open_avg += env.canopy_openness[ilayer] * ca_layer;
		ca_cumm += ca_layer;

	}
	assert(fabs(ca_cumm / G->crown_area - 1) < 1e-6);
	double ca_total = G->crown_area;                   // total crown area
	plant_assim.c_open_avg /= ca_total;                // unitless
	if (by_layer == true){
		plant_assim.dpsi_avg   /= ca_total;                // MPa
		plant_assim.vcmax_avg  /= ca_total;                // umol CO2/m2/s
		plant_assim.gs_avg     /= ca_total;                // mol CO2/m2/s
		plant_assim.vcmax25_avg /= ca_total;               // umol CO2/m2/s
		plant_assim.mc_avg     /= ca_total;                // unitless
		//std::cout << "--- total (by layer) \n";
		//std::cout << "h = " << G->height << ", nz* = " << env.n_layers << ", I = " << plant_assim.c_open_avg << ", fapar = " << fapar << ", A = " << plant_assim.gpp/ca_total << " umol/m2/s x " << ca_total << " = " << plant_assim.gpp << ", vcmax_avg = " << plant_assim.vcmax_avg << "\n"; 
	}

	if (by_layer == false){
		auto res = leaf_assimilation_rate(plant_assim.c_open_avg, fapar, env, par, traits);
		plant_assim.gpp        = (res.a + res.vcmax * par.rd) * ca_total;
		plant_assim.rleaf      = (res.vcmax * par.rd) * ca_total;
		plant_assim.trans      = res.e * ca_total;
		plant_assim.dpsi_avg   = res.dpsi;
		plant_assim.vcmax_avg  = res.vcmax;
		plant_assim.gs_avg     = res.gs;
		plant_assim.vcmax25_avg = res.vcmax25;
		plant_assim.mc_avg     = res.mc;

		// std::cout << "--- total (by avg light)\n";
		// std::cout << "h = " << G->height << ", nz* = " << env.n_layers << ", I = " << plant_assim.c_open_avg << ", fapar = " << fapar << ", A = " << plant_assim.gpp/ca_total << " umol/m2/s x " << ca_total << " = " << plant_assim.gpp << ", vcmax_avg = " << plant_assim.vcmax_avg << "\n"; 
		// std::cout << "a " << res.a << std::endl;
		// std::cout << "vcmax " << res.vcmax << std::endl;
		// std::cout << "rd " << res.rd << std::endl;
		// std::cout << "e " << res.e << std::endl;
		// std::cout << "dpsi " << res.dpsi << std::endl;
		// std::cout << "gs " << res.gs << std::endl;
		// std::cout << "vcmax25 " << res.vcmax25 << std::endl;
		// std::cout << "mc " << res.mc << std::endl;
		// std::cout << "GPP " << plant_assim.gpp << std::endl;
		// std::cout << "rleaf " << plant_assim.rleaf << std::endl;
		// std::cout << "trans " << plant_assim.trans << std::endl;
	}
	// std::cout << "---\nCA traversed = " << ca_cumm << " -- " << G->crown_area << "\n";

	// Convert units from per sec to per unit_t (unit_t is the unit in which time is counted, e.g. yr, day)
	double sec_per_unit_t = 86400 * par.days_per_tunit; // s-1 ---> unit_t-1

	plant_assim.gpp   *= (sec_per_unit_t * 1e-6 * par.cbio);        // umol co2/s ----> umol co2/unit_t --> mol co2/unit_t --> kg/unit_t 
	plant_assim.rleaf *= (sec_per_unit_t * 1e-6 * par.cbio);        // umol co2/s ----> umol co2/unit_t --> mol co2/unit_t --> kg/unit_t 
	plant_assim.trans *= (sec_per_unit_t * 18e-3);                  // mol h2o/s  ----> mol h2o/unit_t  --> kg h2o /unit_t

	// TODO: other traits (vcmax, jmax, gs) etc could also be converted but they are only used in output and not in dynamics
}


template<class Env>
PlantAssimilationResult Assimilator::net_production(Env& env, PlantArchitecture* G, PlantParameters& par, PlantTraits& traits){
	plant_assim = PlantAssimilationResult(); // reset plant_assim

	calc_plant_assimilation_rate(env, G, par, traits); // update plant_assim
	les_update_lifespans(G->lai, par, traits);


	// std::cout << "GPP: " << plant_assim.gpp << std::endl;
	// std::cout << "trans: " << plant_assim.trans << std::endl;
	plant_assim.rleaf = leaf_respiration_rate(G, par, traits);      // kg unit_t-1  
	// std::cout << "Rleaf " << plant_assim.rleaf << std::endl;
	plant_assim.rroot = root_respiration_rate(G, par, traits);     // kg unit_t-1
	// std::cout << "Rroot " << plant_assim.rroot << std::endl;
	plant_assim.rstem = sapwood_respiration_rate(G, par, traits);  // kg unit_t-1
	// std::cout << "Rstem " << plant_assim.rstem << std::endl;
// 
	plant_assim.tleaf = leaf_turnover_rate(kappa_l, G, par, traits);  // kg unit_t-1
	// std::cout << "tleaf " << plant_assim.tleaf << std::endl;
	plant_assim.troot = root_turnover_rate(kappa_r, G, par, traits);  // kg unit_t-1
	// std::cout << "troot " << plant_assim.troot << std::endl;

	//SUPER QUICK FIX!!! 

	// if(std::isnan(plant_assim.gpp)){
	// 	plant_assim.gpp = 0;
	// }
	// if(std::isnan(plant_assim.trans)){
	// 	plant_assim.trans = 0;
	// }
	// if(std::isnan(plant_assim.rroot)){
	// 	plant_assim.rroot = 0;
	// }
	// if(std::isnan(plant_assim.rstem)){	
	// 	plant_assim.rstem = 0;
	// }
	// if(std::isnan(plant_assim.tleaf)){	
	// 	plant_assim.tleaf = 0;
	// }
	// if(std::isnan(plant_assim.troot)){	
	// 	plant_assim.troot = 0;
	// }
	

	double A = plant_assim.gpp;
	double R = plant_assim.rleaf + plant_assim.rroot + plant_assim.rstem;
	double T = plant_assim.tleaf + plant_assim.troot;

	plant_assim.npp = par.y * (A - R) - T; // net biomass growth rate (kg unit_t-1)
	if (std::isnan(plant_assim.npp)){
		std::cout << "ART = " << A << " " << R << " " << T << '\n';
		throw std::runtime_error("assimilator npp is nan");
	}

	// if (G->height > 15) std::cout << "h/A = " << G->height << " / " << A/G->crown_area << std::endl;
	// if (env.n_layers > 1 && G->height < 5) std::cout << "h/L/ml/mr | A/R/T/Vc = " << G->height << " / " << G->lai << " / " << G->leaf_mass(traits) << " / " << G->root_mass(traits) << " | " << A << " / " << R << " / " << T << " / " << plant_assim.vcmax_avg << "\n"; 
	// std::cout.flush();
	return plant_assim;
}

} // namespace plant




