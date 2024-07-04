// cppimport
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <vector>
// #include "pybind_override.cpp"
// #include "environment_base.h"
#include "plantfate_patch.cpp"
#include "light_environment.cpp"
#include "climate.cpp"
#include "pspm_interface.cpp"
#include "plant_architecture.cpp"
#include "assimilation.cpp"
#include "plant.cpp"
#include "state_restore.cpp"
#include "adaptive_species.h"
#include "community_properties.h"
#include "community_properties.cpp"
#include "traits_params.cpp"
#include "climate_stream.cpp"

namespace py = pybind11;
namespace pf = pfate;
namespace pfenv = pfate::env;

PYBIND11_MODULE(pypfate, m)
{
	py::class_<pfenv::Clim>(m, "Clim", py::dynamic_attr())
		.def(py::init<>())
		.def("print", &pfenv::Clim::print)
		.def_readwrite("tc", &pfenv::Clim::tc)
		.def_readwrite("ppfd", &pfenv::Clim::ppfd)
		.def_readwrite("rn", &pfenv::Clim::rn)
		.def_readwrite("vpd", &pfenv::Clim::vpd)
		.def_readwrite("co2", &pfenv::Clim::co2)
		.def_readwrite("elv", &pfenv::Clim::elv)
		.def_readwrite("swp", &pfenv::Clim::swp)
		.def_readwrite("vwind", &pfenv::Clim::vwind)
		.def_readwrite("pa", &pfenv::Clim::pa);

	py::class_<pfenv::Climate>(m, "Climate")
		.def(py::init<>())
		.def("set_elevation", &pfenv::Climate::set_elevation)
		.def("init_co2", &pfenv::Climate::init_co2)
		.def("init_forcing_acclim", &pfenv::Climate::init_forcing_acclim)
		.def("set_acclim_timescale", &pfenv::Climate::set_acclim_timescale)
		.def("set_forcing_acclim", &pfenv::Climate::set_forcing_acclim)
		.def_readwrite("clim_inst", &pfenv::Climate::clim_inst)
		.def_readwrite("clim_acclim", &pfenv::Climate::clim_acclim);

	py::class_<pf::PlantFateConfig>(m, "PlantFateConfig")
		.def(py::init<>())
		.def_readwrite("paramsFile", &pf::PlantFateConfig::paramsFile)
		.def_readwrite("out_dir", &pf::PlantFateConfig::out_dir)
		.def_readwrite("parent_dir", &pf::PlantFateConfig::parent_dir)
		.def_readwrite("expt_dir", &pf::PlantFateConfig::expt_dir)
		.def_readwrite("save_state", &pf::PlantFateConfig::save_state)
		.def_readwrite("state_outfile", &pf::PlantFateConfig::state_outfile)
		.def_readwrite("config_outfile", &pf::PlantFateConfig::config_outfile)
		.def_readwrite("continueFrom_stateFile", &pf::PlantFateConfig::continueFrom_stateFile)
		.def_readwrite("continueFrom_configFile", &pf::PlantFateConfig::continueFrom_configFile)
		.def_readwrite("continuePrevious", &pf::PlantFateConfig::continuePrevious)
		.def_readwrite("saveStateInterval", &pf::PlantFateConfig::saveStateInterval)
		.def_readwrite("n_species", &pf::PlantFateConfig::n_species)
		.def_readwrite("traits_file", &pf::PlantFateConfig::traits_file)
		.def_readwrite("evolve_traits", &pf::PlantFateConfig::evolve_traits)
		.def_readwrite("y0", &pf::PlantFateConfig::y0)
		.def_readwrite("yf", &pf::PlantFateConfig::yf)
		.def_readwrite("ye", &pf::PlantFateConfig::ye)
		.def_readwrite("timestep", &pf::PlantFateConfig::timestep)
		.def_readwrite("T_cohort_insertion", &pf::PlantFateConfig::T_cohort_insertion)
		.def_readwrite("T_seed_rain_avg", &pf::PlantFateConfig::T_seed_rain_avg)
		.def_readwrite("T_return", &pf::PlantFateConfig::T_return)
		.def_readwrite("T_invasion", &pf::PlantFateConfig::T_invasion)
		.def_readwrite("res", &pf::PlantFateConfig::res)
		.def_readwrite("solver_method", &pf::PlantFateConfig::solver_method)
		.def_readwrite("evolvable_traits", &pf::PlantFateConfig::evolvable_traits)
		.def_readwrite("trait_variances", &pf::PlantFateConfig::trait_variances)
		.def_readwrite("trait_scalars", &pf::PlantFateConfig::trait_scalars)
		.def_readwrite("T_r0_avg", &pf::PlantFateConfig::T_r0_avg);

	py::class_<plant::PlantTraits>(m, "PlantTraits")
		.def(py::init<>())
		.def("print", &plant::PlantTraits::print)
		.def_readwrite("lma", &plant::PlantTraits::lma)
		.def_readwrite("zeta", &plant::PlantTraits::zeta)
		.def_readwrite("fcr", &plant::PlantTraits::fcr)
		.def_readwrite("hmat", &plant::PlantTraits::hmat)
		.def_readwrite("fhmat", &plant::PlantTraits::fhmat)
		.def_readwrite("seed_mass", &plant::PlantTraits::seed_mass)
		.def_readwrite("wood_density", &plant::PlantTraits::wood_density)
		.def_readwrite("p50_xylem", &plant::PlantTraits::p50_xylem)
		.def_readwrite("K_leaf", &plant::PlantTraits::K_leaf)
		.def_readwrite("K_xylem", &plant::PlantTraits::K_xylem)
		.def_readwrite("b_leaf", &plant::PlantTraits::b_leaf)
		.def_readwrite("b_xylem", &plant::PlantTraits::b_xylem)
		.def_readwrite("p50_leaf", &plant::PlantTraits::p50_leaf);

	py::class_<plant::PlantParameters>(m, "PlantParameters")
		.def(py::init<>())
		// .def("print", &plant::PlantParameters::print)
		.def_readwrite("cD0", &plant::PlantParameters::cD0)
		.def_readwrite("cD1", &plant::PlantParameters::cD1)
		.def_readwrite("m_alpha", &plant::PlantParameters::m_alpha)
		.def_readwrite("m_beta", &plant::PlantParameters::m_beta)
		.def_readwrite("m_gamma", &plant::PlantParameters::m_gamma);

	py::class_<pf::CommunityProperties> communityProperties(m, "CommunityProperties");

	communityProperties.def(py::init<>())
		.def_readwrite("fluxes", &pf::CommunityProperties::fluxes)
		.def_readwrite("structure", &pf::CommunityProperties::structure)
		.def_readwrite("species", &pf::CommunityProperties::species)
		.def_readwrite("misc", &pf::CommunityProperties::misc)
		.def_readwrite("acc_traits", &pf::CommunityProperties::acc_traits);

	py::class_<pf::CommunityProperties::Fluxes>(communityProperties, "Fluxes")
		.def(py::init<>())
		.def_readwrite("gpp", &pf::CommunityProperties::Fluxes::gpp)
		.def_readwrite("npp", &pf::CommunityProperties::Fluxes::npp)
		.def_readwrite("trans", &pf::CommunityProperties::Fluxes::trans)
		.def_readwrite("gs", &pf::CommunityProperties::Fluxes::gs)
		.def_readwrite("tleaf", &pf::CommunityProperties::Fluxes::tleaf)
		.def_readwrite("troot", &pf::CommunityProperties::Fluxes::troot)
		.def_readwrite("rleaf", &pf::CommunityProperties::Fluxes::rleaf)
		.def_readwrite("rroot", &pf::CommunityProperties::Fluxes::rroot)
		.def_readwrite("rstem", &pf::CommunityProperties::Fluxes::rstem)
		.def_readwrite("mort", &pf::CommunityProperties::Fluxes::mort);

	py::class_<pf::CommunityProperties::Structure>(communityProperties, "Structure")
		.def(py::init<>())
		.def_readwrite("leaf_mass", &pf::CommunityProperties::Structure::leaf_mass)
		.def_readwrite("stem_mass",&pf::CommunityProperties::Structure::stem_mass)
		.def_readwrite("croot_mass", &pf::CommunityProperties::Structure::croot_mass)
		.def_readwrite("froot_mass",&pf::CommunityProperties::Structure::froot_mass)
		.def_readwrite("biomass", &pf::CommunityProperties::Structure::biomass)
		.def_readwrite("basal_area", &pf::CommunityProperties::Structure::basal_area)
		.def_readwrite("canopy_area", &pf::CommunityProperties::Structure::canopy_area)
		.def_readwrite("canopy_area_uc", &pf::CommunityProperties::Structure::canopy_area_uc)
		.def_readwrite("n_ind", &pf::CommunityProperties::Structure::n_ind)
		.def_readwrite("height", &pf::CommunityProperties::Structure::height)
		.def_readwrite("lai", &pf::CommunityProperties::Structure::lai)
		.def_readwrite("lai_vert", &pf::CommunityProperties::Structure::lai_vert);

	py::class_<pf::CommunityProperties::CommunitySpecies>(communityProperties, "CommunitySpecies")
		.def(py::init<>())
		.def_readwrite("n_ind_vec", &pf::CommunityProperties::CommunitySpecies::n_ind_vec)
		.def_readwrite("biomass_vec", &pf::CommunityProperties::CommunitySpecies::biomass_vec)
		.def_readwrite("basal_area_vec", &pf::CommunityProperties::CommunitySpecies::basal_area_vec)
		.def_readwrite("canopy_area_vec", &pf::CommunityProperties::CommunitySpecies::canopy_area_vec)
		.def_readwrite("height_vec", &pf::CommunityProperties::CommunitySpecies::height_vec)
		.def_readwrite("mortality_vec", &pf::CommunityProperties::CommunitySpecies::mortality_vec);

	py::class_<pf::CommunityProperties::Misc>(communityProperties, "Misc")
		.def(py::init<>())
		.def_readwrite("cc_est", &pf::CommunityProperties::Misc::cc_est);

	py::class_<pf::CommunityProperties::Acc_traits>(communityProperties, "")
		.def(py::init<>())
		.def_readwrite("vcmax", &pf::CommunityProperties::Acc_traits::vcmax)
		.def_readwrite("dpsi", &pf::CommunityProperties::Acc_traits::dpsi);

	py::class_<pf::Patch>(m, "Patch")
		.def(py::init<std::string>())
		.def("init", &pf::Patch::init)
		.def("simulate_to", &pf::Patch::simulate_to)
		.def("update_climate", static_cast<void (pf::Patch::*)(double, double, double, double, double)>(&pf::Patch::update_climate))
		.def("update_climate_acclim", &pf::Patch::update_climate_acclim)
		.def("simulate", &pf::Patch::simulate)
		.def("close", &pf::Patch::close)
		.def_readwrite("config", &pf::Patch::config)
		.def_readonly("props", &pf::Patch::props);
};