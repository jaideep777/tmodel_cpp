#ifndef PLANT_FATE_PLANT_ASSIMILATION_H_
#define PLANT_FATE_PLANT_ASSIMILATION_H_

#include <phydro.h>

#include "traits_params.h"
#include "plant_architecture.h"

namespace plant{

/// @brief   Set of all variables calculated by the assimilator
/// @ingroup physiology
struct PlantAssimilationResult{
	double gpp = 0;          ///< Gross plant-level production [kg-biomass unit_t-1] 
	double npp = 0;          ///< Net plant-level production [kg-biomass unit_t-1]
	double trans = 0;        ///< Transpiration [kg-h2o unit_t-1]

	double dpsi_avg = 0;     ///< Soil-leaf water potential difference \f$\Delta\psi\f$
	double vcmax_avg = 0;    ///< Crown-area weighted average Vcmax across canopy layers [umol m-2 s-1] 
	double vcmax25_avg = 0;  ///< Average Vcmax at 25 degC [umol m-2 s-1] 
	double mc_avg = 0;       ///< \f$m_c = (\chi c_a - \Gamma^*)/(\chi c_a + K_M) \f$
	double gs_avg = 0;       ///< Crown-area weighted average stomatal conductance across canopy layers
	double c_open_avg = 0;   ///< Crown-area weighted average canopy opennness experience by the plant

	double rleaf = 0;        ///< Leaf dark respiration rate [kg-biomass unit_t-1]
	double rroot = 0;        ///< Fine root respiration rate [kg-biomass unit_t-1]
	double rstem = 0;        ///< Sapwood respiration rate (excluding coarse root) [kg-biomass unit_t-1]

	double tleaf = 0;        ///< Leaf turnover rate [kg-biomass unit_t-1]
	double troot = 0;        ///< Fine root turnover rate [kg-biomass unit_t-1]
};


/// @brief Object for computing gross and net production, respiration, and turnover.
/// @ingroup physiology
class Assimilator{
	public:
	// ~~ Last calculated phydro values
	// ~~ These are defined here rather than in local scope for debugging purposes. 
	PlantAssimilationResult plant_assim; ///< Plant assimilation result calculated by calc_plant_assimilation_rate()
	// ~~

	double kappa_l;   ///< leaf turnover rate, updated by les functions
	double kappa_r;   ///< fine root turnover rate, updated by les functions

	public:

	/// @brief  Calculate leaf-level assimilation rate using the Phydro model	/// @tparam _Climate 
	/// @param fipar   fraction of top-canopy PAR incident on the crown of this plant
	/// @param fapar   fraction of incident PAR absorbed by the crown
	/// @param clim    forcing variables
	/// @param par     plant params
	/// @param traits  plant traits
	/// @return        leaf assimilatio rate and a bunch of other leaf-level things
	template<class _Climate>
	phydro::PHydroResult leaf_assimilation_rate(double fipar, double fapar, _Climate& clim, PlantParameters& par, PlantTraits& traits);


	/// @brief  Calculate whole-plant gross assimilation, transpiration, gs, etc. 
	template<class Env>
	void  calc_plant_assimilation_rate(Env& env, PlantArchitecture* G, PlantParameters& par, PlantTraits& traits);


	/// @brief  Calculate whole-plant net assimilation 
	template<class Env>
	PlantAssimilationResult net_production(Env& env, PlantArchitecture* G, PlantParameters& par, PlantTraits& traits);


	/// @brief Leaf economics - calculate optimal leaf lifespan 
	/// @{
	void   les_update_lifespans(double lai, PlantParameters& par, PlantTraits& traits);
	double les_assim_reduction_factor(phydro::PHydroResult& res, PlantParameters& par);
	/// @}


	/// @brief Calculate leaf and fine-root respiration rates 
	/// @{
	// leaf respiration rate - should be calculated AFTER asimialtion (needs updated Phydro outputs)
	double leaf_respiration_rate(PlantArchitecture* G, PlantParameters& par, PlantTraits& traits);
	double root_respiration_rate(PlantArchitecture* G, PlantParameters& par, PlantTraits& traits);
	double sapwood_respiration_rate(PlantArchitecture* G, PlantParameters& par, PlantTraits& traits);
	/// @}


	/// @brief Calculate leaf and fine-root turnover rates 
	/// @{
	double leaf_turnover_rate(double _kappa_l, PlantArchitecture* G, PlantParameters& par, PlantTraits& traits);
	double root_turnover_rate(double _kappa_r, PlantArchitecture* G, PlantParameters& par, PlantTraits& traits);
	/// @}

};

} // namespace plant

#include "assimilation.tpp"

#endif


