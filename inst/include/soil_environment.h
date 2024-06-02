#ifndef PLANT_FATE_SOIL_ENVIRONMENT_H_
#define PLANT_FATE_SOIL_ENVIRONMENT_H_

#include <iostream>
#include <vector>

#include <phydro.h>

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

	double swc = 300;    /// soil water content [mm]

	public:
	void spinup(int y, std::vector<double>& sw_in, std::vector<double>& tair, std::vector<double>& pn, std::vector<double>& snow);

	/// @brief Updates soil water content from water fluxes over interval dt
	/// 
	/// @param dt     time interval [days] 
	/// @param precip precipitation [mm day-1]
	/// @param pe_soil potential soil evaporation (i.e., if the soil were saturate with water) [mm day-1]
	void water_balance(double dt, double precip, double pe_soil);

	/// @brief  get soil water potential corresponding to current swc
	/// @return soil water potential [MPa]
	double get_swp();

};

} // env
} // namespace pfate


#endif

