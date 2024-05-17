#ifndef PLANT_FATE_ENV_CLIMATE_STREAM_H_
#define PLANT_FATE_ENV_CLIMATE_STREAM_H_

#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <cmath>
#include <vector>

#include "climate.h"
#include <csvstream.h> // from Flare

namespace pfate{
namespace env{

class ClimateStream{
	public:
	flare::CsvStream i_met_stream;
	flare::CsvStream a_met_stream;
	flare::CsvStream co2_stream;

	std::string i_metFile = "";
	std::string a_metFile = "";
	std::string co2File = "";

	bool update_i_met = false;
	bool update_a_met = false;
	bool update_co2 = false;
	bool use_precip_data = false;

	private:
	template<class T>
	T as(std::string s);

	public:
	void init();
	void updateClimate(double julian_time, Climate& C);
};


template<class T>
T ClimateStream::as(std::string s){
	std::stringstream sin(s);
	T data;
	sin >> data;
	return data;
}


} // namespace env
} // namespace pfate

#endif
