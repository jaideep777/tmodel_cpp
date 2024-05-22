#ifndef PLANT_FATE_PFATE_COMMUNITY_PROPERTIES_H_
#define PLANT_FATE_PFATE_COMMUNITY_PROPERTIES_H_

#include "pspm_interface.h"
#include "adaptive_species.h"
#include "utils/sequence.h"

#ifndef M_PI
#define M_PI 3.14159265358
#endif

namespace pfate{

class Patch;

/// @brief This class calculates and stores community-level properties of interest
class CommunityProperties{
	
	public: 
	// CO2 and water fluxes
	struct Fluxes{
		double gpp = 0;      ///< Gross primary productivity [kgC m-2 day-1]
		double npp = 0;      ///< Net primary productivity [kgC m-2 day-1]
		double trans = 0;    ///< Transpiration [kg m-2 day-1 ~= mm day-1]
		double gs = 0;       ///< Stomatal conductance [mol-h2o m-2 s-1]
		double tleaf = 0;    ///< Leaf turnover rate [kgC m-2 day-1]
		double troot = 0;    ///< Fine root turnover rate [kgC m-2 day-1]
		double rleaf = 0;    ///< Leaf respiration rate [kgC m-2 day-1]
		double rroot = 0;    ///< Fine root respiration rate [kgC m-2 day-1]
		double rstem = 0;    ///< Stem + Coarse-root respiration rate [kgC m-2 day-1]
		double mort = 0;     ///< Biomass loss rate by mortality [kgC m-2 day-1]
	};

	// Structural properties of the community
	struct Structure{
		double leaf_mass = 0;       ///< Leaf biomass [kgC m-2]
		double stem_mass = 0;       ///< Stem biomass [kgC m-2]
		double croot_mass = 0;      ///< Coarse root biomass [kgC m-2]
		double froot_mass = 0;      ///< Fine root biomass [kgC m-2]
		double biomass = 0;         ///< Total biomass [kgC m-2]
		double basal_area = 0;      ///< Total basal area at breast height [m2 m2-1] 
		double canopy_area = 0;     ///< Total crown area [m2 m-2]
		double canopy_area_uc = 0;  ///< Total crown area in the uppermost canopy layer [m2 m-2]
		double n_ind = 0;           ///< Number of individuals in the community [m-2]
		double height = 0;          ///< Average height of all individuals [m]
		double lai = 0;             ///< Community leaf area index [m2 m-2]
		std::vector <double> lai_vert; ///< Vertical profile of LAI (cumulative LAI above height i metres, where i is index in this vector)
	};

	// Species level totals of some properties in the same units as above
	struct CommunitySpecies{
		std::vector<double> n_ind_vec;
		std::vector<double> biomass_vec;
		std::vector<double> basal_area_vec;
		std::vector<double> canopy_area_vec;
		std::vector<double> height_vec;
		std::vector<double> mortality_vec;
		// std::vector<double> mortality_vec_base;
		// std::vector<double> mortality_vec_light;
		// std::vector<double> mortality_vec_hyd;
	};

	// Miscellaneous properties
	struct Misc{
		double cc_est = 0; ///< Estimated carbon cost of leaves averaged over the leaf lifespan
	};

	// Acclimated traits
	struct Acc_traits{
		double vcmax = 0; ///< Average Vcmax in leaves of the upper canopy [umol m-2 s-1]
		double dpsi = 0;  ///< Average dpsi in leaves of the upper canopy [MPa]
	};


	public:
	bool b_output_cohort_props = false;
	std::string emgProps_file;
	std::string cwmAvg_file;
	std::string cwmperSpecies_file;
	std::string traits_file;

	// public:

	Structure structure;
	Fluxes fluxes;
	CommunitySpecies species;
	Misc misc;
	Acc_traits acc_traits;

	void openStreams(std::string dir);
	void closeStreams();
	void writeOut(double t, Patch& P);

	void resize(int n);
	bool isResident(Species_Base* spp);
	void update(double t, Patch& P);

	CommunityProperties& operator /= (double s);
	CommunityProperties& operator += (const CommunityProperties& s);

	private:
	std::vector<std::string> varnames = { "diameter", "height", "lai", "mort", "fec", "rgr", "gpp" };
	std::ofstream cohort_props_out;
	std::ofstream size_dists_out;
	std::ofstream fzst;
	std::ofstream fco;
	std::ofstream flai;
	std::ofstream foutd;
	std::ofstream fouty;
	std::ofstream fouty_spp;
	std::ofstream ftraits;
	std::ofstream fclim;


	template<class Func>
	double integrate_prop(double t, Solver& S, Func f){
		double x = 0;
		for(int k = 0; k < S.n_species(); ++k){
			bool is_resident = static_cast<AdaptiveSpecies<PSPM_Plant>*>(S.species_vec[k])->isResident;
			if(is_resident){
				x += S.state_integral([&S, k, f](int i, double t){
					auto p = (static_cast<Species<PSPM_Plant>*>(S.species_vec[k]))->getCohort(i);
					PSPM_Plant* pp = &p;
					return f(pp);
					}, t, k);
			}
		}
		return x;
	}

	template<class Func>
	std::vector<double> integrate_prop_per_species(double t, Solver& S, Func f){
		std::vector<double> x;
		x.resize(S.n_species(), 0);
		for(int k = 0; k < S.n_species(); ++k){
			bool is_resident = static_cast<AdaptiveSpecies<PSPM_Plant>*>(S.species_vec[k])->isResident;
			if(is_resident){
				x[k] = S.state_integral([&S, k, f](int i, double t){
					auto p = (static_cast<Species<PSPM_Plant>*>(S.species_vec[k]))->getCohort(i);
					PSPM_Plant* pp = &p;
					return f(pp);
					}, t, k);
			}
		}
		return x;
	}

	template<class Func>
	std::vector<double> integrate_prop_above_per_species(double t, double xlow, Solver& S, Func f){
		std::vector<double> x;
		x.resize(S.n_species(), 0);
		for(int k = 0; k < S.n_species(); ++k){
			bool is_resident = static_cast<AdaptiveSpecies<PSPM_Plant>*>(S.species_vec[k])->isResident;
			if(is_resident){
				x[k] = S.integrate_wudx_above([&S, k, f](int i, double t){
					auto p = (static_cast<Species<PSPM_Plant>*>(S.species_vec[k]))->getCohort(i);
					PSPM_Plant* pp = &p;
					return f(pp);
					}, t, { xlow }, k);
			}
		}
		return x;
	}

	template<class Func>
	double integrate_prop_above(double t, double xlow, Solver& S, Func f){
		double x = 0;
		for(int k = 0; k < S.n_species(); ++k){
			bool is_resident = static_cast<AdaptiveSpecies<PSPM_Plant>*>(S.species_vec[k])->isResident;
			if(is_resident){
				x += S.integrate_wudx_above([&S, k, f](int i, double t){
					auto p = (static_cast<Species<PSPM_Plant>*>(S.species_vec[k]))->getCohort(i);
					PSPM_Plant* pp = &p;
					return f(pp);
					}, t, { xlow }, k);
			}
		}
		return x;
	}


};

CommunityProperties operator + (CommunityProperties lhs, CommunityProperties& rhs);

} // namespace pfate

#endif
