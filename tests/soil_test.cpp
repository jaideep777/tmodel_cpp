#include <iostream>
#include <fstream>

#include "climate_stream.h"

using namespace std;

int main(){
	
	env::ClimateStream C;

	C.i_metFile = "tests/data/MetData_AmzFACE_Monthly_2000_2015_PlantFATE_new.csv";
	C.a_metFile = "tests/data/MetData_AmzFACE_Monthly_2000_2015_PlantFATE_new.csv";
	C.co2File = "tests/data/CO2_AMB_AmzFACE2000_2100.csv";
	
	C.update_i_met = true;
	C.update_a_met = true;
	C.update_co2 = true;

	C.init();
	// C.print();

	// C.print_all();

//	C.updateClimate(2001.92);
	env::Climate climate;
	climate.set_acclim_timescale(60);

	ofstream fout("climate_inst.txt");
	ofstream fouta("climate_acclim.txt");
	// for (double t = 1921; t < 2081; t += 1/120.0){
	// for (double t = 2000; t < 2005; t += 1/120.0){
	for (double t = -1000+1e-6; t < 2021+1e-6; t += 1/12.0){
		int year = int(t);
		double month = (t-int(t))*12;
		cout << setprecision(12) << "t = " << t << " id = " << C.i_met_stream.julian_to_indices(flare::yearsCE_to_julian(t)).idx << " (" << year << "/" << month << ")\n";
		C.updateClimate(flare::yearsCE_to_julian(t)+1e-6, climate);
		// C.print_line(t);
		fout << t << "\t" << climate.clim_inst.tc << "\t" << climate.clim_inst.vpd << "\t" << climate.clim_inst.ppfd << "\t" << climate.clim_inst.swp << "\t" << climate.clim_inst.co2 << "\n";
		cout << t << "\t" << climate.clim_inst.tc << "\t" << climate.clim_inst.vpd << "\t" << climate.clim_inst.ppfd << "\t" << climate.clim_inst.swp << "\t" << climate.clim_inst.co2 << "\n";
		fouta << t << "\t" << climate.clim_acclim.tc << "\t" << climate.clim_acclim.vpd << "\t" << climate.clim_acclim.ppfd << "\t" << climate.clim_acclim.swp << "\t" << climate.clim_acclim.co2 << "\n";
	}
	fout.close();

	return 0;

}

//dat = read.delim("~/codes/tmodel_cpp/climate.txt", header=F)
//dato = read.csv("~/codes/tmodel_cpp/tests/data/MetData_AmzFACE_Monthly_2000_2001_PlantFATE.csv", header=T)
//dato$t = dato$Year + (dato$Month-1)/12
//dato = rbind(dato, dato, dato, dato, dato, dato)
//dato$Year = c(rep(2000,12), rep(2001,12), 
//              rep(2002,12), rep(2003,12),
//              rep(2004,12), rep(2005,12),
//              rep(2006,12), rep(2007,12),
//              rep(2008,12), rep(2009,12),
//              rep(2010,12), rep(2011,12))
//dato$t = dato$Year + (dato$Month-1)/12

//par(mfrow=c(4,1), mar=c(4,4,1,1), oma=c(1,1,1,1))
//plot(dat$V2~dat$V1, type="l", col="black")
//points(dato$Temp~dato$t, pch=20, type="p", col="red")
//plot(dat$V3~dat$V1, type="l", col="black")
//points(dato$VPD~dato$t, pch=20, type="p", col="red")
//plot(dat$V4~dat$V1, type="l", col="black")
//points(dato$PAR~dato$t, pch=20, type="p", col="red")
//plot(dat$V5~dat$V1, type="l", col="black")
//points(dato$SWP~dato$t, pch=20, type="p", col="red")
