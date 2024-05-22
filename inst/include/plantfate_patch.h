#ifndef PLANT_FATE_PFATE_PATCH_H_
#define PLANT_FATE_PFATE_PATCH_H_

#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include <numeric>
#include <functional>

#include <solver.h>
#include <time_stepper.h>
#include "pspm_interface.h"
#include "community_properties.h"
#include "adaptive_species.h"
#include "state_restore.h"
#include "climate_stream.h"
#include "plantfate_config.h"

namespace pfate{

class Patch{
	public:
	PlantFateConfig config;

	double t_next_disturbance;
	double t_next_invasion;
	double t_next_savestate;
	double t_next_writestate;

	plant::PlantParameters par0;
	plant::PlantTraits traits0;

	// should be moved out of patch
	env::ClimateStream  climate_stream;
	flare::TimeStepper  ts;

	// io::Initializer     I;
	Solver              S;
	PSPM_Environment    E;

	// SolverIO      sio;
	// SpeciesProps  cwm;
	CommunityProperties props;

	public:
	Patch(std::string params_file);

	void set_i_metFile(std::string file);
	void set_a_metFile(std::string file);
	void set_co2File(std::string co2file);

	void init(double tstart, double tend);

	void simulate_to(double t);

	void update_climate(double julian_time, env::ClimateStream& c_stream);

	void update_climate(double co2, double tc, double vpd, double ppfd, double swp);
	void update_climate_acclim(double t_julian, double co2, double tc, double vpd, double ppfd, double swp);

	void simulate();
	void simulateClimate();

	void close();

	void save(std::ostream& fout);
	void restore(std::istream& fin);

	private:
	double runif(double rmin = 0, double rmax = 1);

	std::vector<plant::PlantTraits> readTraitsFromFile(std::string fname);

	// FIXME: Setting const input seed rain for mutants doesnt work. Is that a problem? 
	/// @brief     Calculate growth rates of all species and update seed input
	/// @param t   Current time
	/// @param dt  timestep over which growth rate is to be calculated
	/// @param S   Solver
	/// @ingroup   trait_evolution
	/// @details   Species growth rate is defined from the seed perspective, i.e.,
	///            \f[r = \frac{1}{\Delta t}log\left(\frac{S_\text{out}}{S_\text{in}}\right),\f] where \f$S\f$ is the seed rain (rate of seed production summed over all individuals of the species)
	void validate_seedrain(std::vector<double>& seeds);
	void calc_r0(double t);

	void removeSpeciesAndProbes(AdaptiveSpecies<PSPM_Plant>* spp);
	void addSpeciesAndProbes(double t, const plant::PlantTraits& traits);
	void shuffleSpecies();

	void removeDeadSpecies(double t);
	void addRandomSpecies(double t);
	void evolveTraits(double t, double dt_evolution);
	void disturbPatch(double t);
};

} // namespace pfate

#endif

