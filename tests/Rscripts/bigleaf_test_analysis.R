library(tidyverse)
rm(list=ls())


input_dir = "~/codes/Plant-FATE/tests/data/"
output_dir = "~/codes/Plant-FATE/pspm_output_test"
# output_dir = "~/output_data/test_3spp_100yr"
expt_dir = "test_2spp_evol_p50" #_old_params"
expt_dir = "test_bigleaf_gfguy" #_old_params"



setwd(paste0(output_dir,"/",expt_dir))


plot_to_file = F


# seeds1 = read.delim("seeds.csv", header=F, col.names = paste0("V", 1:(n_species+2)))
dat = read.csv("bigleaf_D_PFATE.csv", header=T, sep="\t")

dat = dat %>% mutate(date = as.POSIXct(i, orders = "%Y-%m-%d %H:%M:%S"))

if (plot_to_file) png("master_plot.png", width=2412*1.5, height = 1472*1.5, res=300)

par(mfrow=c(4,4), mar=c(4,5,1,1), oma=c(1,1,1,1))

matplot(y=cbind(dat$ppfd_acclim, dat$ppfd_inst), x=dat$date, type="l", lty=1, col=c(scales::muted("seagreen"), "seagreen"), ylab="PPFD\n(umol m-2 s-1)", xlab="Time (years)")
matplot(y=cbind(dat$tc_acclim, dat$tc_inst), x=dat$date, type="l", lty=1, col=c(scales::muted("orange"), "orange"), ylab="Temp\n(degC)", xlab="Time (years)")
matplot(y=cbind(dat$rn_acclim, dat$rn_inst), x=dat$date, type="l", lty=1, col=c(scales::muted("yellow2"), "yellow2"), ylab="Rn\n(W m-2)", xlab="Time (years)")
matplot(y=cbind(dat$vpd_acclim, dat$vpd_inst), x=dat$date, type="l", lty=1, col=c(scales::muted("skyblue"), "skyblue"), ylab="VPD\n(Pa)", xlab="Time (years)")
matplot(y=cbind(dat$swp_acclim, dat$swp_inst), x=dat$date, type="l", lty=1, col=c(scales::muted("cyan3"), "cyan3"), ylab="psi_s\n(Pa)", xlab="Time (years)")
matplot(y=cbind(dat$co2_acclim, dat$co2_inst), x=dat$date, type="l", lty=1, col=c(scales::muted("pink2"), "pink2"), ylab="CO2\n(ppm)", xlab="Time (years)")
matplot(y=cbind(dat$precip_acclim, dat$precip_inst), x=dat$date, type="l", lty=1, col=c(scales::muted("blue"), "blue"), ylab="Precip\n(mm/day)", xlab="Time (years)")
matplot(y=cbind(dat$swc), x=dat$date, type="l", lty=1, col=c("mediumspringgreen"), ylab="Soil WC\n(mm)", xlab="Time (years)")

# par(mfrow=c(4,4), mar=c(4,5,1,1), oma=c(1,1,1,1))

matplot(y=cbind(dat$dpsi), x=dat$date, type="l", lty=1, col=c("seagreen"), ylab="Dpsi\n(MPa)", xlab="Time (years)")
matplot(y=cbind(dat$assim_gross), x=dat$date, type="l", lty=1, col=c("green4", "green3", "brown"), ylab="GPP\n(gC/m2/day)", xlab="Time (years)")
matplot(y=cbind(dat$gs), x=dat$date, type="l", lty=1, col=c("cyan3"), ylab="Stomatal conductance\n(mol/m2/s)", xlab="Time (years)")
matplot(y=cbind(dat$vcmax), x=dat$date, type="l", lty=1, col=c("green3"), ylab="Vcmax\n(umol/m2/s)", xlab="Time (years)")
matplot(y=cbind(dat$jmax), x=dat$date, type="l", lty=1, col=c("green3"), ylab="Jmax\n(umol/m2/s)", xlab="Time (years)")
matplot(y=cbind(dat$transpiration), x=dat$date, type="l", lty=1, col=c("blue"), ylab="T\n(mm/day)", xlab="Time (years)")
matplot(y=cbind(dat$swc), x=dat$date, type="l", lty=1, col=c("mediumspringgreen"), ylab="Soil WC\n(mm)", xlab="Time (years)")
matplot(y=cbind(dat$latent_energy/86400, dat$le_wet_soil/86400), x=dat$date, type="l", lty=1, col=c("green4","yellow3"), ylab="LE\n(J m-2 s-1)", xlab="Time (years)")


if (plot_to_file) dev.off()

# # Convert data to currently usable Plant-FATE format
# 
# dat_inst = read.csv("~/codes/Drought_MIP/drivers_data/PF-CWM_GYF_DY_inst_forcing_daytime_means.csv")
# dat_acclim = read.csv("~/codes/Drought_MIP/drivers_data/PF-CWM_GYF_DY_acclim_forcing_3hr_maxima.csv")
# outfile = "~/codes/Drought_MIP/drivers_data/PF-CWM_GYF_DY_combined_PlantFATE.csv"
# 
# dat_inst %>%
#   select(time, date, sw24hrmean, vpd, temp.C., precip.mm.) %>%
#   rename(sw_day = sw24hrmean) %>%
#   left_join(
#     dat_acclim %>%
#       select(date, shortwave.W.m2.)
#   ) %>%
#   mutate(
#     doy = yday(date),
#   ) %>%
#   reframe(
#     year = year(date),
#     month = month(date),
#     Decimal_year = year + doy/365.2425,
#     temp = temp.C.,
#     vpd = vpd/100,
#     par = sw_day * 2.04,
#     par_max = shortwave.W.m2. * 2.04,
#     swp = 0.05,
#     precip = precip.mm.
#   ) %>%
#   write.csv(file = outfile, row.names = F)



