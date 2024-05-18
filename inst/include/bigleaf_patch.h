#ifndef PLANT_FATE_PFATE_BIGLEAF_PATCH_H_
#define PLANT_FATE_PFATE_BIGLEAF_PATCH_H_

#include <vector>
#include <ostream>

#include "traits_params.h"
#include "assimilation.h"

#include "climate.h"
#include "climate_stream.h"
#include "soil_environment.h"

#include <time_stepper.h>

namespace pfate{


class BigLeafPatch{
	public:
	flare::TimeStepper ts;
	env::ClimateStream climate_stream;

	env::Climate forcing;
	env::SoilEnvironment soil_env;
	double fapar = 1;

	plant::PlantParameters par0;
	plant::PlantTraits traits0;

	io::Initializer I;
	plant::Assimilator bigleaf_assimilator;

	phydro::PHydroResult phydro_out;

	double timestep;
	double t0, tf;

	std::string outfile;

	public:

	BigLeafPatch(std::string params_file);

	void set_i_metFile(std::string file);
	void set_a_metFile(std::string file);
	void set_co2File(std::string co2file);
	void init_co2(double _co2);
	void set_fapar(double _fapar);

	void init(double _t0, double _tf);

	void update_climate(double julian_time);

	std::vector<std::string> get_header();
	void printHeader(std::ostream& lfout);

	std::vector<double> get_state(double t);
	void printState(double t, std::ostream& lfout);

	void printMeta();

	void simulate();
};

} // namespace pfate

#endif
