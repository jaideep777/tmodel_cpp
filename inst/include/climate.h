#ifndef PLANT_FATE_ENV_CLIMATE_H_
#define PLANT_FATE_ENV_CLIMATE_H_

#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <cmath>
#include <vector>
#include <phydro.h>  // for calc_patm()

namespace pfate{
namespace env{

class Clim{
	public:
	double tc = 25.5;         ///< Temperature [C]
	double ppfd = 500;        ///< PAR (daily 24-hr mean) [umol m-2 s-1]
	double rn = 250;          ///< Net radiation at surface [W m-2]
	double vpd = 540;         ///< Vapour pressure deficit [Pa]
	double co2 = 368.9;       ///< Atmospheric CO2 [ppm]
	double elv = 0;           ///< Site elevation [m.a.s.l]
	double swp = -0.04;       ///< Soil water potential [MPa]
	double vwind = 3;         ///< Wind speed [m s-1]
	double pa;                ///< Surface pressure [Pa]
	double precip = 0;        ///< Precipitation [mm day-1]

	Clim();

	void set_elevation(double _elv);

	void print();

	// operators - these are used by exp averager
	Clim& operator += (const Clim& rhs);
	Clim& operator -= (const Clim& s);
	Clim& operator *= (double s);

	private:
	template <class Functor>
	void apply_all(const Clim& rhs, Functor f){
		f(tc, rhs.tc);
		f(ppfd, rhs.ppfd);
		f(rn, rhs.rn);
		f(vpd, rhs.vpd);
		f(co2, rhs.co2);
		f(elv, rhs.elv);
		f(swp, rhs.swp);
		f(vwind, rhs.vwind);
		f(pa, rhs.pa);
		f(precip, rhs.precip);
	}

	template <class Functor>
	void apply_all(double rhs, Functor f){
		f(tc, rhs);
		f(ppfd, rhs);
		f(rn, rhs);
		f(vpd, rhs);
		f(co2, rhs);
		f(elv, rhs);
		f(swp, rhs);
		f(vwind, rhs);
		f(pa, rhs);
		f(precip, rhs);
	}

};

Clim operator + (Clim lhs, const Clim& rhs);
Clim operator - (Clim lhs, const Clim& rhs);
Clim operator * (Clim lhs, double s);
Clim operator * (double s, Clim lhs);


// ClimateStream can probably be merged with this
class Climate{
	protected:
	// These default settings allow the first update to take the new values
	double tau_acclim = 1e-12;   ///< Acclimation timescale [days] 
	double t_last = -1e20;       ///< Time of last updated clim_acclim [julian day] 

	public:
	virtual ~Climate(){};  // need virtual destructor since print() is virtual

	Clim clim_inst;    ///< Mean climate over the timestep
	Clim clim_acclim;  ///< Daily climate measured as mean over 3 hrs around max radiation

	void set_elevation(double _elv);
	void init_co2(double _co2);

	/// @brief Initialize acclim forcing 
	/// @param t0  initial time [julian day]
	/// @param c0  initial forcing value
	void init_forcing_acclim(double t0, const Clim& c0);

	void reset_forcing_acclim();

	/// @brief Set acclimation timescale
	/// @param tau timescale [days]
	void set_acclim_timescale(double tau);

	/// @brief Update acclim forcing as exp-weighted average of current and new value 
	/// @param t0  current time [julian day]
	/// @param c0  current forcing value
	void set_forcing_acclim(double t, const Clim& c);

	virtual void print(double t);
	void print_line(double t);
};

} // namespace env
} // namespace pfate

#endif
