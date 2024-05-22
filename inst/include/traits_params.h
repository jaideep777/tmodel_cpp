#ifndef PLANT_FATE_PLANT_TRAITS_PARAMS_H_
#define PLANT_FATE_PLANT_TRAITS_PARAMS_H_

#include "utils/initializer_v2.h"
#include <io_utils.h>
#include <cmath>
#include <string>
#include <io_utils.h>
#include <cassert>


namespace plant{

/// \ingroup physiology
class PlantTraits{
	protected:
	static std::unordered_map<std::string, double PlantTraits::*> members;

	public:
	std::string species_name = "Tectona grandis";

	// fixed (genetic) traits
	public:
	double lma;             ///< leaf mass per leaf area [kg/m2]
	double zeta;            ///< root mass per leaf area [kg/m2]
	double fcr;             ///< coarse root mass per unit stem mass [-] 
	double hmat;            ///< height at maturity [m]
	double fhmat;           ///< height at reproductive maturity as fraction of hmat
	double seed_mass;       ///< [kg]
	double wood_density;    ///< [kg/m3]
	double p50_xylem;       ///< Xylem hydraulic vulnerability [MPa]
	double K_leaf;          ///< Leaf conductivity [m]
	double K_xylem;         ///< Leaf conductivity [m]
	double b_leaf;          ///< Shape parameter of leaf vulnerabilty curve [-]
	double b_xylem;         ///< Shape parameter of leaf vulnerabilty curve [-]
	double sm_xylem;        ///< Xylem safety margin (P50x - Pg88) [MPa]
	double m;               ///< Crown shape flatness at the top
	double n;               ///< Crown top-heaviness 

	// traits set by coordination
	double a;               ///< Initial height to diameter ratio 
	double c;               ///< Crown area to sapwood area ratio
	double p50_leaf;        ///< Leaf or whole-plant hydraulic vulnerability [MPa] (calculated from Xylem P50 and Safety margin)


	public:
	void init(io::Initializer& I);
	void initFromFile(std::string fname);

	// function to get/set traits by name. Needed for specifying evolvable traits from config file
	double& operator[](const std::string& s);

	// Just for debugging purposes - to check if 2 plants have the same traits
	bool operator == (const PlantTraits& rhs) const;

	void save(std::ostream& fout);
	void restore(std::istream& fin);
	void print();

};


/// \ingroup physiology
class PlantParameters{
	public:
	// Photosynthesis paramaters  
	double kphio;           ///< Quantum use efficiency
	double alpha;           ///< Cost of maintaining photosynthetic capacity
	double gamma;           ///< Cost of hydraulic risks

	// Allocation and geometric paramaters  
	double fg;		        ///< upper canopy gap fraction

	// LAI optimization
	double Cc;                  ///< leaf construction costs
	double Chyd;                ///< hydraulic costs
	double response_intensity;	///< speed of response to environment
	double max_alloc_lai;       ///< max fraction of NPP that can be allocated to LAI increment
	double dl;	                ///< stepsize for profit derivative
	double lai0;                ///< initial lai
	bool   optimize_lai;

	// Leaf Economics
	double les_u;           ///< LES u [dimensionless]
	double les_cc;          ///< Cost of leaf construction and maintenance  [dimensionless]
	double les_k1;          ///< Conversion factor: g biomass / mol CO2 (see cbio below)
	double les_k2;          ///< Conversion factor: (mol-CO2/day) / (umol-CO2/s)
	double les_hT_dH;       ///< Arrhenius temperature response of  [J mol-1]
	double les_hT_c;        ///<   # - 
	double les_molar_R;     ///< Universal gas constant [J mol-1 K-1]

	// Respiration and turnover 
	double rd;              ///< leaf dark respiration rate per unit photosynthetic capacity (r_leaf = rl*vcmax*canopy_area) [yr-1]
	double rr;              ///< fine-root respiration rate [m2 kg-1]
	double rs;              ///< sapwood respiration rate [yr-1]
	double p50x_cost;       ///< Extra respiratory cost of maintining high P50x [yr-1 MPa-2]

	double cbio;            ///< Biomass expansion factor: kg biomass per mol CO2 
	double y;               ///< Growth respiration factor [-]

	double k_light;		    ///< light extincttion coefficient

	// Demographics
	double a_f1;            ///< max fractional allocation to reproduction
	double a_f2;            ///< rate of increase in reproductive investment

	double ll_seed;         ///< longevity of seeds in the seed pool

	// Dispersal and germination
	double Sd;              ///< probability of survival during dispersal
	double npp_Sghalf;      ///< required productivity for 0.5 probability of survival during germination

	// Mortality	
	double cD0, cD1, eD0;
	double m_alpha, m_beta, m_gamma;
	double eWD_alpha, eWD_gamma;
	double cWD0, eWD;
	double m_hydraulic;

	// time stepping
	// below parameter is not explicitly specified in ini file. It is set bythe time stepper for PlantFATE and manually in TreeLife
	double days_per_tunit;      ///< Days in time unit used (e.g., if the unit of time is years, this will be 365.2425)
	double years_per_tunit_avg; ///< Average number of years per time unit, assuming 1 year = 365.2425 days

	public:
	void init(io::Initializer& I);
	void initFromFile(std::string fname);
	void set_tscale(double tscale);
	void print();
	void save(std::ostream& fout);
	void restore(std::istream& fin);

};



} // namespace plant

#endif

