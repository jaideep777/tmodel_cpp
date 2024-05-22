#include "traits_params.h"

namespace plant{

std::unordered_map<std::string, double PlantTraits::*> PlantTraits::members{
	{"lma", &PlantTraits::lma},
	{"zeta", &PlantTraits::zeta},
	{"fcr", &PlantTraits::fcr},
	{"hmat", &PlantTraits::hmat},
	{"fhmat", &PlantTraits::fhmat},
	{"seed_mass", &PlantTraits::seed_mass},
	{"wood_density", &PlantTraits::wood_density},
	{"p50_xylem", &PlantTraits::p50_xylem},
	{"K_leaf", &PlantTraits::K_leaf},
	{"K_xylem", &PlantTraits::K_xylem},
	{"b_leaf", &PlantTraits::b_leaf},
	{"b_xylem", &PlantTraits::b_xylem},
	{"sm_xylem", &PlantTraits::sm_xylem},
	{"m", &PlantTraits::m},
	{"n", &PlantTraits::n},
	{"a", &PlantTraits::a},
	{"c", &PlantTraits::c}
};

void PlantTraits::init(io::Initializer& I){
	lma = I.get<double>("lma");
	zeta = I.get<double>("zeta");
	fcr = I.get<double>("fcr");
	hmat = I.get<double>("hmat");
	fhmat = I.get<double>("fhmat");
	seed_mass = I.get<double>("seed_mass");
	wood_density = I.get<double>("wood_density");
	p50_xylem = I.get<double>("p50_xylem");
	K_leaf = I.get<double>("K_leaf");
	K_xylem = I.get<double>("K_xylem");
	b_leaf = I.get<double>("b_leaf");
	b_xylem = I.get<double>("b_xylem");
	sm_xylem = I.get<double>("sm_xylem");
	m = I.get<double>("m");
	n = I.get<double>("n");
	a = I.get<double>("a");
	c = I.get<double>("c");
	// p50_leaf = // set by coordination
}

void PlantTraits::initFromFile(std::string fname){
	io::Initializer I;
	I.parse(fname);
	init(I);
}

double& PlantTraits::operator[](const std::string& s){
	return this->*(members.find(s)->second);
}

// Just for debugging purposes - to check if 2 plants have the same traits
bool PlantTraits::operator == (const PlantTraits& rhs) const{
	return
		(this->lma == rhs.lma &&
			this->zeta == rhs.zeta &&
			this->fcr == rhs.fcr &&
			this->hmat == rhs.hmat &&
			this->fhmat == rhs.fhmat &&
			this->seed_mass == rhs.seed_mass &&
			this->wood_density == rhs.wood_density &&
			this->p50_xylem == rhs.p50_xylem &&
			this->K_leaf == rhs.K_leaf &&
			this->K_xylem == rhs.K_xylem &&
			this->b_leaf == rhs.b_leaf &&
			this->b_xylem == rhs.b_xylem &&
			this->sm_xylem == rhs.sm_xylem &&
			this->m == rhs.m &&
			this->n == rhs.n &&
			this->a == rhs.a &&
			this->c == rhs.c);
}


// Changelog:
// v2: m,n,a,c moved to traits from parameters
// v3: added p50leaf in save/restore
void PlantTraits::save(std::ostream& fout){
	fout << "Traits::v3   ";
	fout << std::quoted(species_name) << ' ';
	fout << std::make_tuple(
		lma
		, zeta
		, fcr
		, hmat
		, fhmat
		, seed_mass
		, wood_density
		, p50_xylem
		, K_leaf
		, K_xylem
		, b_leaf
		, b_xylem
		, sm_xylem
		, m
		, n
		, a
		, c
		, p50_leaf
	);
	fout << '\n';
}


void PlantTraits::restore(std::istream& fin){
	std::string s; fin >> s; // discard version number
	assert(s == "Traits::v3");

	fin >> std::quoted(species_name);
	fin >> lma
		>> zeta
		>> fcr
		>> hmat
		>> fhmat
		>> seed_mass
		>> wood_density
		>> p50_xylem
		>> K_leaf
		>> K_xylem
		>> b_leaf
		>> b_xylem
		>> sm_xylem
		>> m
		>> n
		>> a
		>> c
		>> p50_leaf;
}

void PlantTraits::print(){
	std::cout << "Traits:\n";
	std::cout << "   lma          = " << lma << '\n';
	std::cout << "   zeta         = " << zeta << '\n';
	std::cout << "   fcr          = " << fcr << '\n';
	std::cout << "   hmat         = " << hmat << '\n';
	std::cout << "   fhmat        = " << fhmat << '\n';
	std::cout << "   seed_mass    = " << seed_mass << '\n';
	std::cout << "   wood_density = " << wood_density << '\n';
	std::cout << "   p50_xylem    = " << p50_xylem << '\n';
	std::cout << "   K_leaf       = " << K_leaf << '\n';
	std::cout << "   K_xylem      = " << K_xylem << '\n';
	std::cout << "   b_leaf       = " << b_leaf << '\n';
	std::cout << "   b_xylem      = " << b_xylem << '\n';
	std::cout << "   sm_xylem     = " << sm_xylem << '\n';
	std::cout << "   m            = " << m << '\n';
	std::cout << "   n            = " << n << '\n';
	std::cout << "   a            = " << a << '\n';
	std::cout << "   c            = " << c << '\n';
	std::cout << "   p50_leaf     = " << p50_leaf << '\n';
}



void PlantParameters::init(io::Initializer& I){
	// #define GET(x) x = I.get<double>(#_x);
	kphio              = I.get<double>("kphio");
	alpha              = I.get<double>("alpha");
	gamma              = I.get<double>("gamma");
	fg                 = I.get<double>("fg");

	Cc                 = I.get<double>("Cc");
	Chyd               = I.get<double>("Chyd");
	response_intensity = I.get<double>("response_intensity");
	max_alloc_lai      = I.get<double>("max_alloc_lai");
	dl                 = I.get<double>("lai_deriv_step");
	lai0               = I.get<double>("lai0");
	optimize_lai       = (I.get<double>("optimize_lai") == 1) ? true : false;

	les_u              = I.get<double>("les_u");
	les_cc             = I.get<double>("les_cc");
	les_k1             = I.get<double>("les_k1");
	les_k2             = I.get<double>("les_k2");
	les_hT_dH          = I.get<double>("les_hT_dH");
	les_molar_R        = I.get<double>("les_molar_R");
	les_hT_c           = I.get<double>("les_hT_c");

	rd                 = I.get<double>("rd");
	rr                 = I.get<double>("rr");
	rs                 = I.get<double>("rs");
	p50x_cost          = I.get<double>("p50x_cost");

	cbio               = I.get<double>("cbio");
	y                  = I.get<double>("y");
	k_light            = I.get<double>("k_light");
	a_f1               = I.get<double>("a_f1");
	a_f2               = I.get<double>("a_f2");

	Sd                 = I.get<double>("Sd");
	npp_Sghalf         = I.get<double>("npp_Sghalf");

	cD0                = I.get<double>("cD0");
	eD0                = I.get<double>("eD0");
	cD1                = I.get<double>("cD1");
	m_alpha            = I.get<double>("m_alpha");
	m_beta             = I.get<double>("m_beta");
	m_gamma            = I.get<double>("m_gamma");
	eWD_alpha          = I.get<double>("eWD_alpha");
	eWD_gamma          = I.get<double>("eWD_gamma");
	cWD0               = I.get<double>("cWD0");
	eWD                = I.get<double>("eWD");
	m_hydraulic        = I.get<double>("m_hydraulic");
}


void PlantParameters::initFromFile(std::string fname){
	io::Initializer I;
	I.parse(fname);
	init(I);
}


void PlantParameters::set_tscale(double tscale){
	days_per_tunit = tscale;
	years_per_tunit_avg = tscale / 365.2425;
}


void PlantParameters::print(){
	std::cout << "Params:\n";
	std::cout
		<< "   kphio = " << kphio << '\n'
		<< "   alpha = " << alpha << '\n'
		<< "   gamma = " << gamma << '\n'
		<< "   fg = " << fg << '\n'
		<< "   Cc = " << Cc << '\n'
		<< "   Chyd = " << Chyd << '\n'
		<< "   response_intensity = " << response_intensity << '\n'
		<< "   max_alloc_lai = " << max_alloc_lai << '\n'
		<< "   dl = " << dl << '\n'
		<< "   lai0 = " << lai0 << '\n'
		<< "   optimize_lai = " << optimize_lai << '\n'
		<< "   les_u = " << les_u << '\n'
		<< "   les_cc = " << les_cc << '\n'
		<< "   les_k1 = " << les_k1 << '\n'
		<< "   les_k2 = " << les_k2 << '\n'
		<< "   les_hT_dH = " << les_hT_dH << '\n'
		<< "   les_molar_R = " << les_molar_R << '\n'
		<< "   les_hT_c = " << les_hT_c << '\n'
		<< "   rd = " << rd << '\n'
		<< "   rr = " << rr << '\n'
		<< "   rs = " << rs << '\n'
		<< "   p50x_cost = " << p50x_cost << '\n'
		<< "   cbio = " << cbio << '\n'
		<< "   y = " << y << '\n'
		<< "   k_light = " << k_light << '\n'
		<< "   a_f1 = " << a_f1 << '\n'
		<< "   a_f2 = " << a_f2 << '\n'
		<< "   ll_seed = " << ll_seed << '\n'
		<< "   Sd = " << Sd << '\n'
		<< "   npp_Sghalf = " << npp_Sghalf << '\n'
		<< "   cD0 = " << cD0 << '\n'
		<< "   eD0 = " << eD0 << '\n'
		<< "   cD1 = " << cD1 << '\n'
		<< "   m_alpha = " << m_alpha << '\n'
		<< "   m_beta = " << m_beta << '\n'
		<< "   m_gamma = " << m_gamma << '\n'
		<< "   eWD_alpha = " << eWD_alpha << '\n'
		<< "   eWD_gamma = " << eWD_gamma << '\n'
		<< "   cWD0 = " << cWD0 << '\n'
		<< "   eWD = " << eWD << '\n'
		<< "   m_hydraulic = " << m_hydraulic << '\n'
		<< "   days_per_tunit = " << days_per_tunit << '\n'
		<< "   years_per_tunit_avg = " << years_per_tunit_avg << '\n';
}

// Changelog:
// v2: m,n,a,c - move to traits from parameters
// v3: add m_hydraulic
// v4: add days_per_t_unit etc
// v5: add p50x_cost
void PlantParameters::save(std::ostream& fout){
	fout << "Params::v5   ";
	fout << std::make_tuple(
		kphio
		, alpha
		, gamma
		, fg
		, Cc
		, Chyd
		, response_intensity
		, max_alloc_lai
		, dl
		, lai0
		, optimize_lai
		, les_u
		, les_cc
		, les_k1
		, les_k2
		, les_hT_dH
		, les_molar_R
		, les_hT_c
		, rd
		, rr
		, rs
		, p50x_cost
		, cbio
		, y
		, k_light
		, a_f1
		, a_f2
		, ll_seed
		, Sd
		, npp_Sghalf
		, cD0
		, eD0
		, cD1
		, m_alpha
		, m_beta
		, m_gamma
		, eWD_alpha
		, eWD_gamma
		, cWD0
		, eWD
		, m_hydraulic
		, days_per_tunit
		, years_per_tunit_avg
	);
	fout << '\n';
}


void PlantParameters::restore(std::istream& fin){
	std::string s; fin >> s; // discard version number
	assert(s == "Params::v5");

	fin >> kphio
		>> alpha
		>> gamma
		>> fg
		>> Cc
		>> Chyd
		>> response_intensity
		>> max_alloc_lai
		>> dl
		>> lai0
		>> optimize_lai
		>> les_u
		>> les_cc
		>> les_k1
		>> les_k2
		>> les_hT_dH
		>> les_molar_R
		>> les_hT_c
		>> rd
		>> rr
		>> rs
		>> p50x_cost
		>> cbio
		>> y
		>> k_light
		>> a_f1
		>> a_f2
		>> ll_seed
		>> Sd
		>> npp_Sghalf
		>> cD0
		>> eD0
		>> cD1
		>> m_alpha
		>> m_beta
		>> m_gamma
		>> eWD_alpha
		>> eWD_gamma
		>> cWD0
		>> eWD
		>> m_hydraulic
		>> days_per_tunit
		>> years_per_tunit_avg
		;
}

} // namespace plant
