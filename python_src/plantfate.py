from datetime import datetime

import numpy as np
import pandas as pd

from pypfate import Patch as patch
from pypfate import Clim


class Model:
    def __init__(self, param_file):
        self.plantFATE_model = patch(str(param_file))
        self.environment = pd.DataFrame(
            columns=[
                "date",
                "tair",
                "ppfd_max",
                "ppfd",
                "vpd",
                "elv",
                "co2",
                "swp",
                "type",
            ]
        )
        self.emergentProps = pd.DataFrame()
        self.speciesProps = pd.DataFrame()

    def runstep(
        self,
        soil_water_potential,
        vapour_pressure_deficit,
        photosynthetic_photon_flux_density,
        temperature,
    ):
        self.plantFATE_model.update_environment(
            temperature,
            photosynthetic_photon_flux_density * 4,
            photosynthetic_photon_flux_density,
            vapour_pressure_deficit * 1000,
            np.nan,
            np.nan,
            soil_water_potential,
        )
        self.plantFATE_model.simulate_step()

        # self.saveEnvironment()
        # self.saveEmergentProps()

        trans = self.plantFATE_model.props.trans / 365
        # return evapotranspiration, soil_specific_depletion_1, soil_specific_depletion_2, soil_specific_depletion_3
        return trans, 0, 0, 0

    def saveEnvironment(self):
        e = pd.DataFrame(
            {
                "date": [
                    self.plantFATE_model.E.tcurrent,
                    self.plantFATE_model.E.tcurrent,
                ],
                "tair": [
                    self.plantFATE_model.E.weightedAveClim.tc,
                    self.plantFATE_model.E.currentClim.tc,
                ],
                "ppfd_max": [
                    self.plantFATE_model.E.weightedAveClim.ppfd_max,
                    self.plantFATE_model.E.currentClim.ppfd_max,
                ],
                "ppfd": [
                    self.plantFATE_model.E.weightedAveClim.ppfd,
                    self.plantFATE_model.E.currentClim.ppfd,
                ],
                "vpd": [
                    self.plantFATE_model.E.weightedAveClim.vpd,
                    self.plantFATE_model.E.currentClim.vpd,
                ],
                "elv": [
                    self.plantFATE_model.E.weightedAveClim.elv,
                    self.plantFATE_model.E.currentClim.elv,
                ],
                "co2": [
                    self.plantFATE_model.E.weightedAveClim.co2,
                    self.plantFATE_model.E.currentClim.co2,
                ],
                "swp": [
                    self.plantFATE_model.E.weightedAveClim.swp,
                    self.plantFATE_model.E.currentClim.swp,
                ],
                "type": ["WeightedAverage", "Instantaneous"],
            }
        )

        self.environment = pd.concat([self.environment, e])

    def saveEmergentProps(self):
        e = pd.DataFrame(
            {
                "date": [self.plantFATE_model.tcurrent],
                "trans": [self.plantFATE_model.props.trans / 365],
                "gs": [self.plantFATE_model.props.gs],
                "gpp": [self.plantFATE_model.props.gpp * 0.5 / 365 * 1000],
                "lai": [self.plantFATE_model.props.lai],
                "npp": [self.plantFATE_model.props.npp * 0.5 / 365 * 1000],
                "cc_est": [self.plantFATE_model.props.cc_est],
                "croot_mass": [self.plantFATE_model.props.croot_mass * 1000 * 0.5],
                "froot_mass": [self.plantFATE_model.props.froot_mass * 1000 * 0.5],
                "lai_vert": [self.plantFATE_model.props.lai_vert],
                "leaf_mass": [self.plantFATE_model.props.leaf_mass * 1000 * 0.5],
                "resp_auto": [self.plantFATE_model.props.resp_auto * 0.5 / 365 * 1000],
                "stem_mass": [self.plantFATE_model.props.stem_mass * 1000 * 0.5],
            }
        )
        self.emergentProps = pd.concat([self.emergentProps, e])

    def exportEnvironment(self, out_file):
        self.environment.to_csv(out_file, sep=",", index=False, encoding="utf-8")

    def exportEmergentProps(self, out_file):
        self.emergentProps.to_csv(out_file, sep=",", index=False, encoding="utf-8")

    def exportSpeciesProps(self, out_file):
        self.speciesProps.to_csv(out_file, sep=",", index=False, encoding="utf-8")

    def first_step(
        self,
        tstart,
        tend,
        temperature0,  # degrees Celcius, mean temperature
        photosynthetic_photon_flux_density0,
        vapour_pressure_deficit0,
        soil_water_potential0
    ):

        newclim = Clim()
        newclim.tc = temperature0  # C
        newclim.ppfd_max = photosynthetic_photon_flux_density0 * 4
        newclim.ppfd = photosynthetic_photon_flux_density0
        newclim.vpd = vapour_pressure_deficit0 * 1000  # kPa -> Pa
        newclim.swp = soil_water_potential0  # MPa

        self.plantFATE_model.init(tstart, tend, newclim)

    def simulate(self):
        self.plantFATE_model.simulate()

    def calculate_soil_water_potential_MPa(
        self,
        soil_moisture,  # [m]
        soil_moisture_wilting_point,  # [m]
        soil_moisture_field_capacity,  # [m]
        soil_tickness,  # [m]
        wilting_point=-1500,  # kPa
        field_capacity=-33,  # kPa
    ):
        # https://doi.org/10.1016/B978-0-12-374460-9.00007-X (eq. 7.16)
        soil_moisture_fraction = soil_moisture / soil_tickness
        assert soil_moisture_fraction >= 0 and soil_moisture_fraction <= 1
        del soil_moisture
        soil_moisture_wilting_point_fraction = (
            soil_moisture_wilting_point / soil_tickness
        )
        assert (
            soil_moisture_wilting_point_fraction >= 0
            and soil_moisture_wilting_point_fraction <= 1
        )
        del soil_moisture_wilting_point
        soil_moisture_field_capacity_fraction = (
            soil_moisture_field_capacity / soil_tickness
        )
        assert (
            soil_moisture_field_capacity_fraction >= 0
            and soil_moisture_field_capacity_fraction <= 1
        )
        del soil_moisture_field_capacity

        n_potential = -(
            np.log(wilting_point / field_capacity)
            / np.log(
                soil_moisture_wilting_point_fraction
                / soil_moisture_field_capacity_fraction
            )
        )
        assert n_potential >= 0
        a_potential = (
            1.5 * 10**6 * soil_moisture_wilting_point_fraction**n_potential
        )
        assert a_potential >= 0
        soil_water_potential = -a_potential * soil_moisture_fraction ** (-n_potential)
        return soil_water_potential / 1_000_000  # Pa to MPa

    def calculate_vapour_pressure_deficit_kPa(self, temperature, relative_humidity):
        assert (
            temperature < 100
        )  # temperature is in Celsius. So on earth should be well below 100.
        assert (
            temperature > -100
        )  # temperature is in Celsius. So on earth should be well above -100.
        assert (
            relative_humidity >= 1 and relative_humidity <= 100
        )  # below 1 is so rare that it shouldn't be there at the resolutions of current climate models, and this catches errors with relative_humidity as a ratio [0-1].
        # https://soilwater.github.io/pynotes-agriscience/notebooks/vapor_pressure_deficit.html
        saturated_vapour_pressure = 0.611 * np.exp(
            (17.502 * temperature) / (temperature + 240.97)
        )  # kPa
        actual_vapour_pressure = (
            saturated_vapour_pressure * relative_humidity / 100
        )  # kPa
        vapour_pressure_deficit = saturated_vapour_pressure - actual_vapour_pressure
        return vapour_pressure_deficit

    def calculate_photosynthetic_photon_flux_density(self, shortwave_radiation, xi=0.5):
        # https://search.r-project.org/CRAN/refmans/bigleaf/html/Rg.to.PPFD.html
        photosynthetically_active_radiation = shortwave_radiation * xi
        photosynthetic_photon_flux_density = (
            photosynthetically_active_radiation * 4.6
        )  #  W/m2 -> umol/m2/s
        return photosynthetic_photon_flux_density

    def get_plantFATE_input(
        self,
        soil_moisture_layer_1,  # m
        soil_moisture_layer_2,  # m
        soil_moisture_layer_3,  # m
        soil_tickness_layer_1,  # m
        soil_tickness_layer_2,  # m
        soil_tickness_layer_3,  # m
        soil_moisture_wilting_point_1,  # m
        soil_moisture_wilting_point_2,  # m
        soil_moisture_wilting_point_3,  # m
        soil_moisture_field_capacity_1,  # m
        soil_moisture_field_capacity_2,  # m
        soil_moisture_field_capacity_3,  # m
        temperature,  # degrees Celcius, mean temperature
        relative_humidity,  # percentage [0-100]
        shortwave_radiation,  # W/m2, daily mean
    ):
        assert (
            temperature < 100
        )  # temperature is in Celsius. So on earth should be well below 100.
        assert relative_humidity >= 0 and relative_humidity <= 100

        soil_water_potential = self.calculate_soil_water_potential_MPa(
            soil_moisture_layer_1 + soil_moisture_layer_2 + soil_moisture_layer_3,
            soil_moisture_wilting_point_1
            + soil_moisture_wilting_point_2
            + soil_moisture_wilting_point_3,
            soil_moisture_field_capacity_1
            + soil_moisture_field_capacity_2
            + soil_moisture_field_capacity_3,
            soil_tickness_layer_1 + soil_tickness_layer_2 + soil_tickness_layer_3,
        )

        vapour_pressure_deficit = self.calculate_vapour_pressure_deficit_kPa(
            temperature, relative_humidity
        )

        photosynthetic_photon_flux_density = (
            self.calculate_photosynthetic_photon_flux_density(shortwave_radiation)
        )

        return (
            soil_water_potential,
            vapour_pressure_deficit,
            photosynthetic_photon_flux_density,
            temperature,
        )

    def step(
        self,
        soil_moisture_layer_1,  # ratio [0-1]
        soil_moisture_layer_2,  # ratio [0-1]
        soil_moisture_layer_3,  # ratio [0-1]
        soil_tickness_layer_1,  # m
        soil_tickness_layer_2,  # m
        soil_tickness_layer_3,  # m
        soil_moisture_wilting_point_1,  # ratio [0-1]
        soil_moisture_wilting_point_2,  # ratio [0-1]
        soil_moisture_wilting_point_3,  # ratio [0-1]
        soil_moisture_field_capacity_1,  # ratio [0-1]
        soil_moisture_field_capacity_2,  # ratio [0-1]
        soil_moisture_field_capacity_3,  # ratio [0-1]
        temperature,  # degrees Celcius, mean temperature
        relative_humidity,  # percentage [0-100]
        shortwave_radiation,  # W/m2, daily mean
    ):
        (
            soil_water_potential,
            vapour_pressure_deficit,
            photosynthetic_photon_flux_density,
            temperature,
        ) = self.get_plantFATE_input(
            soil_moisture_layer_1,  # ratio [0-1]
            soil_moisture_layer_2,  # ratio [0-1]
            soil_moisture_layer_3,  # ratio [0-1]
            soil_tickness_layer_1,  # m
            soil_tickness_layer_2,  # m
            soil_tickness_layer_3,  # m
            soil_moisture_wilting_point_1,  # ratio [0-1]
            soil_moisture_wilting_point_2,  # ratio [0-1]
            soil_moisture_wilting_point_3,  # ratio [0-1]
            soil_moisture_field_capacity_1,  # ratio [0-1]
            soil_moisture_field_capacity_2,  # ratio [0-1]
            soil_moisture_field_capacity_3,  # ratio [0-1]
            temperature,  # degrees Celcius, mean temperature
            relative_humidity,  # percentage [0-100]
            shortwave_radiation,  # W/m2, daily mean
        )

        (
            evapotranspiration,
            soil_specific_depletion_1,
            soil_specific_depletion_2,
            soil_specific_depletion_3,
        ) = self.runstep(
            soil_water_potential,
            vapour_pressure_deficit,
            photosynthetic_photon_flux_density,
            temperature,
        )

        soil_specific_depletion_1 = (
            np.nan
        )  # this is currently not calculated in plantFATE, so just setting to np.nan to avoid confusion
        soil_specific_depletion_2 = (
            np.nan
        )  # this is currently not calculated in plantFATE, so just setting to np.nan to avoid confusion
        soil_specific_depletion_3 = (
            np.nan
        )  # this is currently not calculated in plantFATE, so just setting to np.nan to avoid confusion

        evapotranspiration = evapotranspiration / 1000  # kg H2O/m2/day to m/day

        return (
            evapotranspiration,
            soil_specific_depletion_1,
            soil_specific_depletion_2,
            soil_specific_depletion_3,
        )

    def finalize(self):
        self.plantFATE_model.close()

    @property
    def n_individuals(self):
        return sum(self.plantFATE_model.cwm.n_ind_vec)

    @property
    def biomass(self):
        return sum(self.plantFATE_model.cwm.biomass_vec)  # kgC / m2
