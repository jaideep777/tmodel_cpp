#include "community_properties.h"
#include "light_environment.h"
#include "plantfate_patch.h"
#include <phydro.h>

using namespace std;

namespace pfate{

void CommunityProperties::resize(int n){
	species.n_ind_vec.resize(n);
	species.biomass_vec.resize(n);
	species.basal_area_vec.resize(n);
	species.canopy_area_vec.resize(n);
	species.height_vec.resize(n);
}

bool CommunityProperties::isResident(Species_Base* spp){
	return static_cast<AdaptiveSpecies<PSPM_Plant>*>(spp)->isResident;
}

void CommunityProperties::update(double t, Patch& P){

	Solver& S = P.S;

	// multiplier to convert unit_t-1 --> day-1
	double m1 = 1 / P.par0.days_per_tunit;

	// multiplier to convert kg biomass --> kg C
	double m2 = 1 / 2.04;

	fluxes.gpp    = m1 * m2 * integrate_prop(t, S, [](PSPM_Plant* p){return p->res.gpp;});
	fluxes.npp    = m1 * m2 * integrate_prop(t, S, [](PSPM_Plant* p){return p->res.npp;});
	fluxes.trans  = m1 * integrate_prop(t, S, [](PSPM_Plant* p){return p->res.trans;});
	fluxes.tleaf  = m1 * m2 * integrate_prop(t, S, [](PSPM_Plant* p){return p->res.tleaf;});
	fluxes.troot  = m1 * m2 * integrate_prop(t, S, [](PSPM_Plant* p){return p->res.troot;});
	fluxes.rleaf  = m1 * m2 * integrate_prop(t, S, [](PSPM_Plant* p){return p->res.rleaf;});
	fluxes.rroot  = m1 * m2 * integrate_prop(t, S, [](PSPM_Plant* p){return p->res.rroot;});
	fluxes.rstem  = m1 * m2 * integrate_prop(t, S, [](PSPM_Plant* p){return p->res.rstem;});
	fluxes.mort   = m1 * m2 * integrate_prop(t, S, [](PSPM_Plant* p){return p->get_biomass() * p->rates.dmort_dt;});
	fluxes.gs = (fluxes.trans * 55.55 / 86400) / 1.6 / (static_cast<PSPM_Environment*>(S.env)->clim_inst.vpd / 1.0325e5); // FIXME: replace hardcoded pressure with clim_inst.pa
	//     ^ convert kg/m2/day --> mol/m2/s

	// Note: for these vector calcs, multipliers are multiplied inside lambdas, where necesary
	species.n_ind_vec       = integrate_prop_above_per_species(t, 0.1, S, [](PSPM_Plant* p){return 1;});
	species.biomass_vec     = integrate_prop_above_per_species(t, 0.1, S, [](PSPM_Plant* p){return p->get_biomass();});
	species.basal_area_vec  = integrate_prop_above_per_species(t, 0.1, S, [](PSPM_Plant* p){
		double D = p->geometry.diameter_at_height(1.3, p->traits);
		return M_PI * D * D / 4;
		});
	species.canopy_area_vec = integrate_prop_above_per_species(t, 0.1, S, [](PSPM_Plant* p){return p->geometry.crown_area;});
	species.height_vec      = integrate_prop_above_per_species(t, 0.1, S, [](PSPM_Plant* p){return p->geometry.height;});
	for (int k=0; k < S.n_species(); ++k) species.height_vec[k] /= species.n_ind_vec[k];
	species.mortality_vec   = integrate_prop_per_species(t, S, [m1, m2](PSPM_Plant* p){return m1 * m2 * p->get_biomass() * p->rates.dmort_dt;});

	structure.n_ind         = std::accumulate(species.n_ind_vec.begin(), species.n_ind_vec.end(), 0.0);
	structure.biomass       = std::accumulate(species.biomass_vec.begin(), species.biomass_vec.end(), 0.0);
	structure.basal_area    = std::accumulate(species.basal_area_vec.begin(), species.basal_area_vec.end(), 0.0);
	structure.canopy_area   = std::accumulate(species.canopy_area_vec.begin(), species.canopy_area_vec.end(), 0.0);

	structure.lai           =      integrate_prop(t, S, [](PSPM_Plant* p){return p->geometry.crown_area * p->geometry.lai;});
	structure.leaf_mass     = m2 * integrate_prop(t, S, [](PSPM_Plant* p){return p->geometry.leaf_mass(p->traits);});
	structure.stem_mass     = m2 * integrate_prop(t, S, [](PSPM_Plant* p){return p->geometry.stem_mass(p->traits);});
	structure.croot_mass    = m2 * integrate_prop(t, S, [](PSPM_Plant* p){return p->geometry.coarse_root_mass(p->traits);});
	structure.froot_mass    = m2 * integrate_prop(t, S, [](PSPM_Plant* p){return p->geometry.root_mass(p->traits);});

	// upper canopy acclimated traits
	double h_uc = static_cast<PSPM_Environment*>(S.env)->z_star[0];
	structure.canopy_area_uc = integrate_prop_above(t, 0.1, S, [h_uc](PSPM_Plant* p){return (p->geometry.height > h_uc) ? p->geometry.crown_area : 0;});

	acc_traits.vcmax = integrate_prop_above(t, 0.1, S, [h_uc](PSPM_Plant* p){return (p->geometry.height > h_uc) ? (p->res.vcmax25_avg * p->geometry.crown_area) : 0; });
	acc_traits.vcmax /= structure.canopy_area_uc;

	acc_traits.dpsi  = integrate_prop_above(t, 0.1, S, [h_uc](PSPM_Plant* p){return (p->geometry.height > h_uc) ? (p->res.dpsi_avg * p->geometry.crown_area) : 0; });
	acc_traits.dpsi /= structure.canopy_area_uc;

	misc.cc_est = (fluxes.tleaf + fluxes.troot + fluxes.rleaf + fluxes.rroot + fluxes.rstem) / fluxes.tleaf;

	// LAI vertical profile
	structure.lai_vert.clear();
	structure.lai_vert.resize(25, 0);
	for (int iz=0; iz < 25; ++iz){
		structure.lai_vert[iz] = integrate_prop(t, S, [iz](PSPM_Plant* p){return p->geometry.crown_area_above(iz, p->traits) * p->geometry.lai;});
	}

	// soil evaporation 
	// TODO: Might be possible to avoid duplication by splitting phydro::env computations
	auto cc = static_cast<PSPM_Environment*>(S.env)->clim_inst;
	double gamma = phydro::calc_psychro(cc.tc, cc.pa);
	double epsilon = phydro::calc_sat_slope(cc.tc) / gamma;
	double lv = phydro::calc_enthalpy_vap(cc.tc);
	double fapar = 1-exp(-0.5*structure.lai);
	fluxes.pe_soil = (1-fapar)*cc.rn*(epsilon/(1+epsilon)) / lv * 86400; // J m-2 s-1 / (J kg-1) * (s/day) = kg m-2 s-1 * s day-1 = kg m-2 day-1
}


CommunityProperties operator + (CommunityProperties lhs, CommunityProperties& rhs){
	lhs += rhs;
	return lhs;
}


void CommunityProperties::openStreams(std::string dir){

	if (b_output_cohort_props){
		cohort_props_out.open(dir + "/cohort_props.csv");
		cohort_props_out << "t,speciesID,cohortID,";
		for (auto vname : varnames) cohort_props_out << vname << ",";
		cohort_props_out << std::endl;
		cohort_props_out << std::setprecision(12);
	}

	size_dists_out.open(dir + "/size_distributions.csv");
	fzst.open(std::string(dir + "/z_star.csv").c_str());
	fco.open(std::string(dir + "/canopy_openness.csv").c_str());
	flai.open(std::string(dir + "/lai_profile.csv").c_str());
	foutd.open(std::string(dir + "/" + emgProps_file).c_str());
	fouty.open(std::string(dir + "/" + cwmAvg_file).c_str());
	fouty_spp.open(std::string(dir + "/" + cwmperSpecies_file).c_str());
	ftraits.open(std::string(dir + "/" + traits_file).c_str());
	// fclim.open(std::string(dir + "/climate_co2.csv").c_str());

	foutd << "YEAR,GPP,NPP,RAU,MORT,GS,ET,PESOIL,VCMAX,DPSI,CCEST,CO2\n";
	fouty << "YEAR,DE,CL,CW,CCR,CFR,CR,CA,BA,TB,LAI\n";
	fouty_spp << "YEAR,PID,DE,PH,CA,BA,TB,SEEDS\n";
	ftraits << "YEAR,SPP,RES,LMA,WD,HMAT,P50X,ZETA,r0_last,r0_avg,r0_exp,r0_cesaro\n";
	// fclim << "t,tc,ppfd_max,ppfd,vpd,co2,elv,swp\n";

}

void CommunityProperties::closeStreams(){
	if (b_output_cohort_props) cohort_props_out.close();
	size_dists_out.close();

	fzst.close();
	fco.close();
	flai.close();
	foutd.close();
	fouty.close();
	fouty_spp.close();
	ftraits.close();
	// fclim.close();
}

void CommunityProperties::writeOut(double t, Patch& P){
	Solver* S = &P.S;

	// multiplier to convert unit_t-1 --> day-1
	double m1 = 1 / P.par0.days_per_tunit;

	// consistently output date in decimal years across all files
	double date = flare::julian_to_yearsCE(P.ts.to_julian(t));

	for (int s=0; s < S->species_vec.size(); ++s){
		auto spp = static_cast<AdaptiveSpecies<PSPM_Plant>*>(S->species_vec[s]);

		std::vector<double> breaks = my_log_seq(0.01, 10, 100);
		std::vector<double> dist = S->getDensitySpecies1D(s, 0, breaks);
		//cout << "here: " << breaks.size() << " " << dist.size() << endl;

		if (spp->isResident){
			size_dists_out << date << "," << spp->species_name << ",";
			for (int i=0; i < 100; ++i){
				size_dists_out << dist[i] << ",";
			}
			size_dists_out << "\n";
		}

		if (b_output_cohort_props){
			if (spp->isResident){
				for (int j=0; j < spp->xsize() - 1; ++j){
					auto& C = spp->getCohort(j);
					cohort_props_out
						<< date << ","
						<< spp->species_name << ","  // use name instead of index s becuase it is unique and order-insensitive
						<< j << ","
						<< C.geometry.diameter << ","
						<< C.geometry.height << ","
						<< C.geometry.lai << ","
						<< C.rates.dmort_dt << ","
						<< C.rates.dseeds_dt << ","
						<< C.rates.rgr << ","
						<< C.res.gpp / C.geometry.crown_area << ",";
					cohort_props_out << "\n";
				}
			}
		}
	}

	foutd
		<< date << ","
		<< fluxes.gpp << ","
		<< fluxes.npp << ","
		<< (fluxes.rleaf + fluxes.rroot + fluxes.rstem) << ","  // kgC/m2/d
		<< fluxes.mort << ","
		<< fluxes.gs << ","
		<< fluxes.trans << ","
		<< fluxes.pe_soil << ","
		<< acc_traits.vcmax << ","
		<< acc_traits.dpsi << ","
		<< misc.cc_est << ","
		<< static_cast<PSPM_Environment*>(S->env)->clim_inst.co2
		<< std::endl;

	fouty
		<< date << ","
		<< structure.n_ind << ","
		<< structure.leaf_mass << ","
		<< structure.stem_mass << ","
		<< structure.croot_mass << ","
		<< structure.froot_mass << ","
		<< (structure.croot_mass + structure.froot_mass) << ","
		<< structure.canopy_area << ","
		<< structure.basal_area << ","
		<< structure.biomass << ","
		<< structure.lai
		<< std::endl;

	for (int k=0; k < S->species_vec.size(); ++k){
		auto spp = static_cast<AdaptiveSpecies<PSPM_Plant>*>(S->species_vec[k]);
		fouty_spp
			<< date << ","
			<< spp->species_name << "," // use name instead of index k becuase it is unique and order-insensitive
			<< species.n_ind_vec[k] << ","
			<< species.height_vec[k] << ","
			<< species.canopy_area_vec[k] << ","
			<< species.basal_area_vec[k] << ","
			<< species.biomass_vec[k] << ","
			<< spp->seeds_hist1.get() * m1
			<< std::endl;
	}

	for (int k=0; k < S->species_vec.size(); ++k){
		auto spp = static_cast<AdaptiveSpecies<PSPM_Plant>*>(S->species_vec[k]);
		ftraits
			<< date << ","
			<< spp->species_name << "," // use name instead of index k becuase it is unique and order-insensitive
			<< spp->isResident << ","
			<< spp->getCohort(-1).traits.lma << ","
			<< spp->getCohort(-1).traits.wood_density << ","
			<< spp->getCohort(-1).traits.hmat << ","
			<< spp->getCohort(-1).traits.p50_xylem << ","
			<< spp->getCohort(-1).traits.zeta << ","
			<< spp->r0_hist.get_last() << ","
			<< spp->r0_hist.get() << ","
			<< 0 << ","
			<< 0 //spp->r0_hist.get_cesaro();
			<< std::endl;
	}

	ftraits.flush();
	fouty_spp.flush();

	flai << date << ",";
	for (int i=0; i < structure.lai_vert.size(); ++i) flai << structure.lai_vert[i] << ",";
	flai << std::endl;

	fzst << date << ",";
	for (auto z : static_cast<PSPM_Environment*>(S->env)->z_star) fzst << z << ",";
	fzst << std::endl;

	fco << date << ",";
	for (auto z : static_cast<PSPM_Environment*>(S->env)->canopy_openness) fco << z << ",";
	fco << std::endl;

}

// FIXME: These operators are not tested
CommunityProperties& CommunityProperties::operator /= (double s){

	fluxes.gpp            /= s;
	fluxes.npp            /= s;
	fluxes.trans          /= s;
	fluxes.gs             /= s;
	fluxes.tleaf          /= s;
	fluxes.troot          /= s;
	fluxes.rleaf          /= s;
	fluxes.rroot          /= s;
	fluxes.rstem          /= s;

	structure.leaf_mass   /= s;
	structure.stem_mass   /= s;
	structure.croot_mass  /= s;
	structure.froot_mass  /= s;
	structure.biomass     /= s;
	structure.basal_area  /= s;
	structure.canopy_area /= s;
	structure.n_ind       /= s;
	structure.height      /= s;
	structure.lai         /= s;
	transform(structure.lai_vert.begin(), structure.lai_vert.end(), structure.lai_vert.begin(), [s](const double& c){ return c / s; });

	misc.cc_est           /= s;

	acc_traits.vcmax      /= s;
	acc_traits.dpsi       /= s;

	transform(species.n_ind_vec.begin(), species.n_ind_vec.end(), species.n_ind_vec.begin(), [s](const double& c){ return c / s; });
	transform(species.biomass_vec.begin(), species.biomass_vec.end(), species.biomass_vec.begin(), [s](const double& c){ return c / s; });
	transform(species.basal_area_vec.begin(), species.basal_area_vec.end(), species.basal_area_vec.begin(), [s](const double& c){ return c / s; });
	transform(species.canopy_area_vec.begin(), species.canopy_area_vec.end(), species.canopy_area_vec.begin(), [s](const double& c){ return c / s; });
	transform(species.height_vec.begin(), species.height_vec.end(), species.height_vec.begin(), [s](const double& c){ return c / s; });

	return *this;
}

CommunityProperties& CommunityProperties::operator += (const CommunityProperties& rhs){

	fluxes.gpp            += rhs.fluxes.gpp;
	fluxes.npp            += rhs.fluxes.npp;
	fluxes.trans          += rhs.fluxes.trans;
	fluxes.gs             += rhs.fluxes.gs;
	fluxes.tleaf          += rhs.fluxes.tleaf;
	fluxes.troot          += rhs.fluxes.troot;
	fluxes.rleaf          += rhs.fluxes.rleaf;
	fluxes.rroot          += rhs.fluxes.rroot;
	fluxes.rstem          += rhs.fluxes.rstem;

	structure.leaf_mass   += rhs.structure.leaf_mass;
	structure.stem_mass   += rhs.structure.stem_mass;
	structure.croot_mass  += rhs.structure.croot_mass;
	structure.froot_mass  += rhs.structure.froot_mass;
	structure.biomass     += rhs.structure.biomass;
	structure.basal_area  += rhs.structure.basal_area;
	structure.canopy_area += rhs.structure.canopy_area;
	structure.n_ind       += rhs.structure.n_ind;
	structure.height      += rhs.structure.height;
	structure.lai         += rhs.structure.lai;
	transform(structure.lai_vert.begin(), structure.lai_vert.end(), rhs.structure.lai_vert.begin(), structure.lai_vert.begin(), std::plus<double>());

	misc.cc_est           += rhs.misc.cc_est;

	acc_traits.vcmax      += rhs.acc_traits.vcmax;
	acc_traits.dpsi       += rhs.acc_traits.dpsi;

	transform(species.n_ind_vec.begin(), species.n_ind_vec.end(), rhs.species.n_ind_vec.begin(), species.n_ind_vec.begin(), std::plus<double>());
	transform(species.biomass_vec.begin(), species.biomass_vec.end(), rhs.species.biomass_vec.begin(), species.biomass_vec.begin(), std::plus<double>());
	transform(species.basal_area_vec.begin(), species.basal_area_vec.end(), rhs.species.basal_area_vec.begin(), species.basal_area_vec.begin(), std::plus<double>());
	transform(species.canopy_area_vec.begin(), species.canopy_area_vec.end(), rhs.species.canopy_area_vec.begin(), species.canopy_area_vec.begin(), std::plus<double>());
	transform(species.height_vec.begin(), species.height_vec.end(), rhs.species.height_vec.begin(), species.height_vec.begin(), std::plus<double>());

	return *this;
}



} // namespace pfate