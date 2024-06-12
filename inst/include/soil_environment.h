#ifndef PLANT_FATE_SOIL_ENVIRONMENT_H_
#define PLANT_FATE_SOIL_ENVIRONMENT_H_

#include <iostream>
#include <vector>

#include <phydro.h>
#include <SPLASH.h>

namespace pfate{
namespace env{

class SoilEnvironment {
	public:
	struct SoilPars{
		double b = 3;        /// Shape parameter of water potential curve [-]
		double rzwsc = 300;  /// Root zone water storage capacity [mm]
		double sstar = 100;  /// global minimum rzwsc [mm]
	} par; /// Soil parameters (may be eventually merged with simulation parameters)

	struct SplashPars{
		double slope;
		double asp;
		double lat;
		double elev;
		std::vector<double> soil_info;
	} par_spl;

	struct SoilState{
		double wn;    // soil moisture [mm]
		double swe;   // snow water equivalent [mm]
		double qin;   // drainage
		double td;    // days left to drain the area upslope
		double nd;    // snow age [days]
	} state;

	double swc = 300;    /// soil water content [mm]

	bool need_spinup = true;

	SPLASH splash;
	SOLAR solar;
	smr dsoil;
	etr dvap;

	public:
	SoilEnvironment();

	void init();

	void spinup(int y, std::vector<double>& sw_in, std::vector<double>& tair, std::vector<double>& pn, std::vector<double>& snow);

	void update_radiation(int doy, int y, double sw_in, double tair, double pn, double snow);

	template <class ForcingInput>
	void update_forcings(int doy, int year, ForcingInput& forcing);

	void water_balance_splash(int doy, int y, double sw_in, double tair, double pn, double snow, double plant_uptake);

	/// @brief Updates soil water content from water fluxes over interval dt
	/// 
	/// @param dt     time interval [days] 
	/// @param precip precipitation [mm day-1]
	/// @param pe_soil potential soil evaporation (i.e., if the soil were saturate with water) [mm day-1]
	void water_balance(double dt, double precip, double pe_soil);

	/// @brief  get soil water potential corresponding to current swc
	/// @return soil water potential [MPa]
	double get_swp();

	void save(std::ostream& fout);
	void restore(std::istream& fin);
};

} // env
} // namespace pfate

#include <soil_environment.tpp>

#endif

