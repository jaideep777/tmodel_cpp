#include <iomanip>
#include <fstream>
#include <random>
#include <string>
#include "plant_params.h"
#include "plant_geometry.h"
#include "assimilation.h"
#include "plant.h"

#include "climate.h"
#include "light_environment.h"

using namespace std;

class Environment : public env::Climate, public env::LightEnvironment {
public:
	void print(double t){
		Climate::print(t);
		LightEnvironment::print();
	}
};

int main(){
	
	// Creating random number generator for soil water potential
	ofstream f_mean("fitness_v_swpmean.txt");
	
	// SET SOIL WATER POTENTIAL SERIES PARAMETERS
	double swp_mean_start = -0.1;
	double swp_mean_stop = -10.0
	double swp_stddev = -4.0;
	double time_step_swp = 10;
	
	// SET SIMULATION PARAMETERS
	double sim_start_t = 2000;
	double sim_end_t = 2100;
	double sim_dt = 0.1;
	
	std::vector<double> swp_mean;
	
	f_mean << "SWP_Mean" << "\t"
	       << "Fitness_Avg" << "\n";
	
	
	for (double i_mean = swp_mean_start; i_mean >= swp_mean_stop; i_mean = i_mean - 0.1 ){
		
		std::vector<double> germi_avg;
		cout << "Starting SWP mean : " << to_string(i_mean) << "\n";
		
		
		// STARTING SIMULATION FOR 100 CYCLES FOR AVERAGING
		for (double iter = 1; iter <=100; iter=iter+1 ){
			
			double germinated = 0;
			
			plant::Plant P;
			P.initParamsFromFile("tests/params/p.ini");
			P.geometry.set_lai(1);
			P.set_size(0.01);
			
			Environment C;
			C.metFile = "tests/data/MetData_AmzFACE_Monthly_2000_2015_PlantFATE.csv";
			C.co2File = "tests/data/CO2_AMB_AmzFACE2000_2100.csv";
			
			std::default_random_engine generator;
			std::normal_distribution<double> dist(i_mean, swp_stddev);
			for (double t=sim_start_t; t<=sim_end_t; t=t+time_step_swp){
				C.t_met.push_back(t);
				double val = dist(generator);
				if(val>0) val = 0;
				C.v_met_swp.push_back(val);
			}
			//C.print(0);
			
			double total_prod = P.get_biomass();
			double total_rep = 0;
			double litter_pool = 0;
			
			for (double t=sim_start_t; t<=sim_end_t; t=t+sim_dt){
				P.grow_for_dt(t, sim_dt, C, total_prod, total_rep, litter_pool, germinated);
			}
			
			germi_avg.push_back(germinated);
			
			if (abs((P.get_biomass()+litter_pool+total_rep)/total_prod - 1) > 2e-5) return 1;
		}
	
		double average = accumulate(germi_avg.begin(), germi_avg.end(), 0.0) / germi_avg.size();
		
		f_mean << i_mean << "\t"
		<< average << "\n";
		germi_avg.clear();
	}
	f_mean.close();
	return 0;
}

