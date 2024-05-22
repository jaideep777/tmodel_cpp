#ifndef UTILS_MATH_MOVING_AVG_H_
#define UTILS_MATH_MOVING_AVG_H_

#include <iostream>

#include <list>
#include <cmath>
#include <cassert>
#include <vector>
#include <algorithm>

class MovingAverager{
	private:
	double area_sum = 0;
	std::list<double> areas;
	std::list<double> dts;
	std::list<double> f_hist;
	std::list<double> t_hist;

	public:
	// seed output history
	double T = 5;
	bool debug = false;

	public:
	
	inline void set_interval(double _T){
		T = _T;
	}
	
	inline void push(double t, double f){
		if (f_hist.empty()){
			f_hist.push_back(f);
			t_hist.push_back(t);
			//areas.push_back(0);
			return;
		}
		
		assert(t > t_hist.back());

		double f_lo = f_hist.back();
		double t_lo = t_hist.back();
		double area_new = (f + f_lo) / 2 * (t - t_lo);
		f_hist.push_back(f);
		t_hist.push_back(t);
		dts.push_back(t - t_lo);
		areas.push_back(area_new);
		area_sum += area_new;
		
		while(!t_hist.empty()){
			if (t_hist.front() < t - T){  // past t is beyond the averaging interval
				area_sum -= areas.front();
				f_hist.pop_front();
				t_hist.pop_front();
				dts.pop_front();
				areas.pop_front();
			}
			else{
				break;
			}
		}
	}

	inline std::size_t size(){
		return t_hist.size();
	}

	inline double get(){
		if (t_hist.size() == 0) return 0;
		else if (t_hist.size() == 1) return f_hist.front();
		else return area_sum/(t_hist.back()-t_hist.front());
	}

	inline double get_first(){
		if (t_hist.size() == 0) return 0;
		else return f_hist.front();
	}

	inline double get_last(){
		if (t_hist.size() == 0) return 0;
		else return f_hist.back();
	}

	inline double get_last_t(){
		if (t_hist.size() == 0) return -1e20;
		else return t_hist.back();
	}

	inline double get_interval(){
		if (t_hist.size() == 0) return 0;
		else return (t_hist.back()-t_hist.front());
	}

	inline void clear(){
		area_sum = 0;
		f_hist.clear();
		t_hist.clear();
		areas.clear();
		dts.clear();
	}
	
	inline void print(){
		std::cout << "MovingAverager:\n";
		auto f_it = f_hist.begin();
		auto t_it = t_hist.begin();
		auto a_it = areas.begin();
		auto d_it = dts.begin();
		
		std::cout << "   t\tf\tA\tdt\n";
		while(f_it != f_hist.end()){
			std::cout << "   " << *t_it << "\t" << *f_it << "\t";
			std::cout << ((a_it != areas.end())? *a_it : 0) << "\t";
			std::cout << ((d_it != dts.end())? *d_it : 0) << "\n";
			++f_it; ++t_it; ++a_it;
		}
		std::cout << "-----\n";
		std::cout << "sum = " << area_sum << ", avg = " << get() << "\n";
		
	}

	inline void print_summary(){
		std::cout << "MovingAverager:  t = " << t_hist.front() << " - " << t_hist.back() << ", x = " << f_hist.front() << " - " << f_hist.back() << ", xmean = " << get() << ", (T = " << t_hist.back() - t_hist.front() << "), npoints = " << f_hist.size() << "\n";
	}

	inline double get_exp(double c = 0){
		if (t_hist.size() == 0) return 0;
		if (t_hist.size() == 1) return f_hist.front();

		if (debug) std::cout << "Exp avg: ";
		auto f_it = f_hist.rbegin();
		auto t_it = t_hist.rbegin();
		double avg = 0;
		double D = 0;
		double t0 = *t_hist.rbegin();
		int i=0;
		while (f_it != prev(f_hist.rend())){
			double t_hi = *t_it;
			double t_lo = *next(t_it);
			double f_hi = (*f_it);
			double w_hi = exp(-c*(t0-t_hi));
			double f_lo = (*next(f_it));
			double w_lo = exp(-c*(t0-t_lo));
			double dt = t_hi - t_lo;
			avg += (f_hi*w_hi + f_lo*w_lo)/2 * dt;
			D += (w_hi + w_lo)/2 * dt;
			++f_it; ++t_it;

			if (debug) std::cout << (f_hi*w_hi + f_lo*w_lo)/2 * dt << " ";
		}
		avg /= D;

		if (debug) std::cout << " | " << avg << std::endl;

		return avg;
	}


	inline double get_cesaro(double c=0){
		if (t_hist.size() == 0) return 0;
		if (t_hist.size() == 1) return f_hist.front();

		// f_hist has 2 or more elements, so compute
		// assert that t is equally spaced
		auto tt = std::minmax_element(dts.begin(), dts.end());
		double thresh = 0.01;
		if (fabs(*tt.first - *tt.second) > thresh) throw std::runtime_error("Cesaro Average requires equally spaced intervals.");
		double Dt = *tt.first;

		// 
		int nmax = f_hist.size();
		std::vector<double> avgs(nmax, 0);

		for (int n=1; n<nmax; ++n){
			avgs[n] = 0;
			double D = 0;
			auto f_it = f_hist.rbegin();
			for (int i=0; i<n; ++i, ++f_it){
				avgs[n] += (*f_it) * exp(-c*i*Dt);
				D += exp(-c*i*Dt);
			}
			avgs[n] /= D;
		}

		double avg = 0;
		double D = 0;
		for (int i=1; i<nmax; ++i){
			avg += avgs[i]*exp(-c*(i-1)*Dt); // Note avgs[0] is not defined, the first number is at index 1
			D += exp(-c*(i-1)*Dt);
		}
		avg /= D;

		if (debug){
			std::cout << "Cesaro: ";
			for (auto x : avgs) std::cout << x << " ";
			std::cout << " | " << avg << std::endl;
		}

		return avg;
	}

	void save(std::ostream &fout){
		fout << "MovingAverager::v1\n";

		fout << area_sum << ' '
		     << T << '\n';
		
		fout << areas.size() << " | "; 
		for (auto x : areas) fout << x << ' '; 
		fout << '\n';

		fout << dts.size() << " | "; 
		for (auto x : dts) fout << x << ' '; 
		fout << '\n';

		fout << f_hist.size() << " | "; 
		for (auto x : f_hist) fout << x << ' '; 
		fout << '\n';

		fout << t_hist.size() << " | "; 
		for (auto x : t_hist) fout << x << ' '; 
		fout << '\n';
		
	}

	void restore(std::istream &fin){
		std::string s; fin >> s; // discard version number
		assert(s == "MovingAverager::v1");

		fin >> area_sum
		    >> T;
		
		int n;
		double x;

		fin >> n >> s; 
		for (int i=0; i<n; ++i) {fin >> x; areas.push_back(x);}

		fin >> n >> s; 
		for (int i=0; i<n; ++i) {fin >> x; dts.push_back(x);}

		fin >> n >> s; 
		for (int i=0; i<n; ++i) {fin >> x; f_hist.push_back(x);}

		fin >> n >> s; 
		for (int i=0; i<n; ++i) {fin >> x; t_hist.push_back(x);}
	}

};


#endif
