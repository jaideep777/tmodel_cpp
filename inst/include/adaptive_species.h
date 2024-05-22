#ifndef PLANT_FATE_PFATE_ADAPTIVE_SPECIES_H_
#define PLANT_FATE_PFATE_ADAPTIVE_SPECIES_H_

#include <solver.h>
#include <vector>
#include <string>
#include <io_utils.h>
#include "utils/moving_average.h"

// Extend the Species class from libpspm to allow trait evolution

/// @defgroup trait_evolution Trait Evolution
/// @brief   This is a collection of all classes and variables that are used in the 
///          evolutionary algorithm. 
///
///          Any changes to the algorithm or to the evolving traits will 
///          typically require changes in these variables.

namespace pfate{

/// @ingroup trait_evolution
/// @brief   This classes extends the Species class provided by the solver
///          to implement the evolutionary algorithm.
/// @tparam  Model 'Model' represents an individual plant which is used as a
///          prototype to construct the indiviuals in the species.
template <class Model>
class AdaptiveSpecies : public Species<Model>{
	private:
	double fg_dx = 0.001;

	public:
	std::string species_name;
	bool isResident = true;
	double t_introduction = 0;    ///< Time since introduction [years]

	double invasion_fitness;
	double r0;
	std::vector<AdaptiveSpecies<Model>*> probes; // /*NO_SAVE_RESTORE*/ This vector is excluded from save/restore because it contains pointers. Species-probes associations are saved/restored by name in Patch::save/restore
	std::vector<std::string> evolvable_traits;
	std::vector<double> fitness_gradient;
	std::vector<double> trait_variance;
	std::vector<double> trait_scalars;     // these scalars will be applied to fg_dx

	MovingAverager seeds_hist1;
	MovingAverager seeds_hist2;
	MovingAverager r0_hist;


	public:
	AdaptiveSpecies(Model M, bool res = true);

	void set_traits(std::vector<double> tvec);
	std::vector<double> get_traits();

	// Species(vector<double> tvec, double u0, bool res);

	void createVariants(Model M);

	void calcFitnessGradient();
	void evolveTraits(double dt);

	void print_extra();

	void save(std::ostream& fout) override;
	void restore(std::istream& fin) override;

};

} // namespace pfate

#include "adaptive_species.tpp"

#endif

