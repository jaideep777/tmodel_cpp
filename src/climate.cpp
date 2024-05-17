#include <iomanip>
#include <cmath>
#include <stdexcept>

#include "climate.h"
using namespace std;

namespace pfate{
namespace env{

Clim::Clim(){
	pa = phydro::calc_patm(elv);
}

void Clim::set_elevation(double _elv){
	elv = _elv;
	pa = phydro::calc_patm(elv);
}


void Clim::print(){
	std::cout << "Clim:\n";
	std::cout << "   tc = " << tc << '\n';
	std::cout << "   ppfd = " << ppfd << '\n';
	std::cout << "   rn = " << rn << '\n';
	std::cout << "   vpd = " << vpd << '\n';
	std::cout << "   co2 = " << co2 << '\n';
	std::cout << "   elv = " << elv << '\n';
	std::cout << "   swp = " << swp << '\n';
	std::cout << "   vwind = " << vwind << '\n';
	std::cout << "   pa = " << pa << '\n';
	std::cout << "   precip = " << precip << '\n';
}

// operators - these are used by exp averager
Clim& Clim::operator += (const Clim& rhs){
	apply_all(rhs, [](double& l, double r){l += r;});
	return *this;
}

Clim& Clim::operator -= (const Clim& s){
	apply_all(s, [](double& l, double r){l -= r;});
	return *this;
}

Clim& Clim::operator *= (double s){
	apply_all(s, [](double& l, double r){l *= r;});
	return *this;
}

Clim operator + (Clim lhs, const Clim& rhs){
	lhs += rhs;
	return lhs;
}

Clim operator - (Clim lhs, const Clim& rhs){
	lhs -= rhs;
	return lhs;
}

Clim operator * (Clim lhs, double s){
	lhs *= s;
	return lhs;
}

Clim operator * (double s, Clim lhs){
	lhs *= s;
	return lhs;
}


void Climate::set_elevation(double _elv){
	clim_inst.set_elevation(_elv);
	clim_acclim.set_elevation(_elv);
}

void Climate::init_co2(double _co2){
	clim_inst.co2 = _co2;
	clim_acclim.co2 = _co2;
}

void Climate::init_forcing_acclim(double t0, const Clim& c0){
	t_last = t0;
	clim_acclim = c0;
}

void Climate::set_acclim_timescale(double tau){
	tau_acclim = tau;
}

// REF: https://stackoverflow.com/questions/1023860/exponential-moving-average-sampled-at-varying-times
void Climate::set_forcing_acclim(double t, const Clim& c){
	double dt = t - t_last;
	double alpha = 1 - exp(-dt / tau_acclim);
	clim_acclim += alpha * (c - clim_acclim);   // Compute exp-weighted moving average
	t_last = t;
}

void Climate::print(double t){
	cout << "Current climate (inst / acclim):\n";
	cout << "   tc      = " << clim_inst.tc << " / " << clim_acclim.tc << '\n';
	cout << "   ppfd    = " << clim_inst.ppfd << " / " << clim_acclim.ppfd << '\n';
	cout << "   vpd     = " << clim_inst.vpd << " / " << clim_acclim.vpd << '\n';
	cout << "   co2     = " << clim_inst.co2 << " / " << clim_acclim.co2 << '\n';
	cout << "   elv     = " << clim_inst.elv << " / " << clim_acclim.elv << '\n';
	cout << "   swp     = " << clim_inst.swp << " / " << clim_acclim.swp << '\n';
	cout << "   pa      = " << clim_inst.pa << " / " << clim_acclim.pa << '\n';
	cout << "   rn      = " << clim_inst.rn << " / " << clim_acclim.rn << '\n';
	cout << "   vwind   = " << clim_inst.vwind << " / " << clim_acclim.vwind << '\n';
	cout << "   precip  = " << clim_inst.precip << " / " << clim_acclim.precip << '\n';
}

void Climate::print_line(double t){
	cout << "Current climate (inst)   (t=" << t << "): ";
	cout << clim_inst.tc << " ";
	cout << clim_inst.ppfd << " ";
	cout << clim_inst.rn << " ";
	cout << clim_inst.vpd << " ";
	cout << clim_inst.co2 << " ";
	cout << clim_inst.elv << " ";
	cout << clim_inst.swp << " ";
	cout << clim_inst.precip << " ";
	cout << '\n';
	cout << "Current climate (acclim) (t=" << t << "): ";
	cout << clim_acclim.tc << " ";
	cout << clim_acclim.ppfd << " ";
	cout << clim_acclim.rn << " ";
	cout << clim_acclim.vpd << " ";
	cout << clim_acclim.co2 << " ";
	cout << clim_acclim.elv << " ";
	cout << clim_acclim.swp << " ";
	cout << clim_acclim.precip << " ";
	cout << '\n';

}

} // namespace env
} // namespace pfate

