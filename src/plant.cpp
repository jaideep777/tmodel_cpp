#include "plant.h"
using namespace std;

namespace plant{

/// TODO: Should init() functions be made constructors, to prevent risk of creating a plant without
/// @brief  This function initializes the plant (traits, par, and geometry) from an Initialzer object
void Plant::init(const PlantParameters& _par, const PlantTraits& _traits){
	par = _par;
	traits = _traits;
	coordinateTraits();
}


/// @brief  This function initializes the plant (traits, par, and geometry) from an Initialzer object
void Plant::init(io::Initializer& I){
	par.init(I);
	traits.init(I);
	init(par, traits);
}

void Plant::initFromFile(std::string file){
	io::Initializer I;
	I.parse(file);
	init(I);
}


void Plant::coordinateTraits(){
	traits.coordinate();
	geometry.init(par, traits);

	// std::cout << "Coordinated traits:\n";
	// traits.print();
}


void Plant::set_size(double x){
	geometry.set_size(x, traits);
}

double Plant::get_biomass() const{
	return geometry.total_mass(traits);
}

void Plant::set_evolvableTraits(const std::vector<std::string>& tnames, const std::vector<double>& tvalues){
	for (int i=0; i < tnames.size(); ++i){
		traits[tnames[i]] = tvalues[i];
	}
	init(par, traits);
}

std::vector<double> Plant::get_evolvableTraits(const std::vector<std::string>& tnames){
	std::vector<double> tvec(tnames.size());
	for (int i=0; i < tnames.size(); ++i){
		tvec[i] = traits[tnames[i]];
	}
	return tvec;
}


void Plant::print(){
	cout << "Plant size:\n";
	cout << "  height = " << geometry.height << "\n";
	cout << "  diameter = " << geometry.diameter << "\n";
	cout << "  crown_area = " << geometry.crown_area << "\n";
	traits.print();
}


}	// namespace plant




