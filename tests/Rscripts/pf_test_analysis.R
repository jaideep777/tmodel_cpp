library(tidyverse)
rm(list=ls())


input_dir = "~/codes/Plant-FATE/tests/data/"
output_dir = "~/codes/Plant-FATE/pspm_output_test"
# output_dir = "~/output_data/test_3spp_100yr"
expt_dir = "test_2spp_evol_p50" #_old_params"
expt_dir = "test_3spp_100yr" #_old_params"

# expt_dir = "test_1spp_evol_p50" #_old_params"
# expt_dir = "test_4spp_evol_wd_hmat" #_old_params"
# expt_dir = "cont_test_main" #_old_params"
solver = "IEBT"

setwd(paste0(output_dir,"/",expt_dir))


plot_to_file = F
plot_trait_space = F

add_band = function(){
  polygon(x=c(2000,5000,5000,2000), y=c(-1e20,-1e20,1e20,1e20), border = NA, col=scales::alpha("yellow2",0.2))
}

add_hband = function(ylim, col=scales::alpha("grey30",0.2), xlim=c(-1e20,3000)){
  polygon(y=c(ylim[1],ylim[2],ylim[2],ylim[1]), x=c(xlim[1],xlim[1],xlim[2],xlim[2]), border = NA, col=col)
}

# seeds1 = read.delim("seeds.csv", header=F, col.names = paste0("V", 1:(n_species+2)))
Zp = read.csv("z_star.csv", header=F, col.names = paste0("V", 1:50))
# BA1 = read.csv("basal_area.csv", header=F, col.names = paste0("V", 1:(n_species+2)))
co = read.csv("canopy_openness.csv", header=F, col.names = paste0("V", 1:50))
lai_v = read.csv("lai_profile.csv", header=F, col.names = paste0("V", 1:27))
traits = read.csv("traits.csv")
dat = read.csv("D_PFATE.csv")
# dat$YEAR = decimal_date(as_date(dat$YEAR, format = "%Y-%m-%d %H:%M:%S GMT (doy = %j)"))
dat2 = read.csv("Y_PFATE.csv")
dat3 = read.csv("Y_mean_PFATE.csv", )
dist = read.csv("size_distributions.csv", header=F)
dist = dist[,-ncol(dist)]
x = exp(seq(log(0.01), log(10), length.out=100))
traits_obs = read.csv(file = paste0(input_dir, "/Amz_trait_filled_HD.csv"))
traits_used = read.csv(file = paste0(input_dir, "/Traits_random_HD2.csv"))

# To get avg size distribution, sum over species and average over years
names(dist)[1:2] = c("YEAR", "SPP")
dist_amb = dist %>% filter(YEAR>min(1101,max(YEAR)-2) & YEAR<2000) %>% filter(YEAR > min(YEAR)) %>% pivot_longer(cols=-(YEAR:SPP), names_to="size_class") %>% group_by(YEAR,size_class) %>% summarize(de = sum(value, na.rm=T)) %>% pivot_wider(names_from = size_class, values_from = de) %>% colMeans(na.rm=T)
dist_ele = dist %>% filter(YEAR>min(4101,max(YEAR)-2) & YEAR<5000) %>% filter(YEAR > min(YEAR)) %>% pivot_longer(cols=-(YEAR:SPP), names_to="size_class") %>% group_by(YEAR,size_class) %>% summarize(de = sum(value, na.rm=T)) %>% pivot_wider(names_from = size_class, values_from = de) %>% colMeans(na.rm=T)

# dist_amb = dist %>% filter(YEAR == 1100) %>% pivot_longer(cols=-(YEAR:SPP), names_to="size_class") %>% group_by(YEAR,size_class) %>% summarize(de = sum(value, na.rm=T)) %>% pivot_wider(names_from = size_class, values_from = de) %>% colMeans(na.rm=T)
# dist_ele = dist %>% filter(YEAR == 1101) %>% pivot_longer(cols=-(YEAR:SPP), names_to="size_class") %>% group_by(YEAR,size_class) %>% summarize(de = sum(value, na.rm=T)) %>% pivot_wider(names_from = size_class, values_from = de) %>% colMeans(na.rm=T)

n_species = dat2 %>% filter(!grepl("probe", .$PID)) %>% pull(PID) %>% unique() %>% length()
n_year = length(unique(dat2$YEAR))
col_species = rainbow(n = n_species, start = 0, end = 0.7, alpha = min(10/n_species, 1))

if (plot_to_file) png("master_plot.png", width=2412*1.5, height = 1472*1.5, res=300)

par(mfcol=c(4,5), mar=c(4.5,6,.5,1), oma=c(1,1,2,1), cex.lab=1.1, cex.axis=1.1, mgp=c(3.2,1,0), las=1)
seeds = dat2 %>% filter(!grepl("probe", .$PID)) %>% select(YEAR, PID, SEEDS) %>% spread(value = "SEEDS", key = "PID")
seeds_smooth = seeds %>% pivot_longer(-YEAR) %>% drop_na() %>% group_by(name) %>% mutate(value = loess(value~YEAR, span=60/length(value)) %>% fitted()) %>% pivot_wider(names_from=name)
matplot(seeds$YEAR, seeds[,-1], lty=1, col=scales::alpha(col_species, 0.3), type="l",
        las=1, xlab="Time (years)", ylab="Species seed\noutput", log="")
# matplot(seeds_smooth$YEAR, seeds_smooth[,-1], lty=1, type="l",
#         col=col_species, #col=scales::alpha(scales::muted(col_species), alpha=min(30/n_species, 1)), 
#         las=1, xlab="Time (years)", ylab="Species seed\noutput", log="")
mtext(line=0.5, side=3, text=expt_dir)
add_band()

# matplot(seeds1$V1, seeds1[,-1], lty=1, col=rainbow(n = n_species+1, start = 0, end = 0.85), type="l",
#         las=1, xlab="Time (years)", ylab="Species Seed output", log="")
# mtext(line=0.5, side=3, text=expt_dir)

BA = dat2 %>% filter(!grepl("probe", .$PID)) %>% select(YEAR, PID, BA) %>% spread(value = "BA", key = "PID")
matplot(BA$YEAR, cbind(BA[,-1], BA %>% select(-YEAR) %>% rowSums(na.rm=T))*1e4, lty=1, col=c(col_species, "black"), type="l",
        las=1, xlab="Time (years)", ylab="Basal area", log="")
add_hband(c(31.29, 31.29*1.02))
add_band()

matplot(Zp$V1, Zp[,-1], lty=1, col=rainbow(n = 10, start = 0, end = 0.85), type="l",
        las=1, xlab="Time (years)", ylab="Z*")
add_band()
# matplot(co$V1, co[,-1], lty=1, col=rainbow(n = 10, start = 0, end = 0.85), type="l",
#         las=1, xlab="Time (years)", ylab="Io")

matplot(y=cbind(dat$DPSI), x=dat$YEAR, type="l", lty=1, col=c("seagreen"), ylab="Dpsi\n(MPa)", xlab="Time (years)")
matlines(y=cbind(fitted(loess(dat$DPSI~dat$YEAR, span=60/n_year))), x=dat$YEAR, type="l", lty=1, col="black", lwd=c(1,0.5))
# add_hband(c(0.16, 0.16555))#, col=scales::alpha("cyan4", 0.6))
# abline(h=c(0.16), col=scales::muted("cyan3"))
add_band()


matplot(y=1:25, x=t(lai_v[,2:26]), lty=1, col=rainbow(n = n_year, start = 0, end = 0.85, alpha=0.05), type="l",
        las=1, xlab="Cumulative LAI", ylab="Height")


plot(dat3$LAI~dat3$YEAR, type="l", col="red3", ylim=c(0,max(dat$LAI,6.5)), xlab="Time (years)", ylab="Total LAI")
# abline(h=c(5.3, 6.2), col=scales::muted("red"))
add_hband(c(5.3, 6.2))#, col=scales::alpha("red3", 0.2))
add_band()
# abline(h=c(3.5), col=scales::muted("grey100"))


matplot(y=cbind(dat$GPP, dat$NPP, dat$MORT)*365.2425, x=dat$YEAR, type="l", lty=1, col=c("green4", "green3", "brown"), ylab="GPP, NPP, MORT\n(kgC/m2/yr)", xlab="Time (years)")
matlines(y=cbind(fitted(loess(dat$GPP~dat$YEAR, span=60/n_year)), 
                 fitted(loess(dat$NPP~dat$YEAR, span=60/n_year)))*365, x=dat$YEAR, type="l", lty=1, col="black", lwd=c(1,0.5))
# points(y=dat$NPP/dat$GPP*4, x=dat$YEAR, type="l", lty=1, col=c("yellow1"))
# abline(h=c(3,3.5), col="grey")
add_hband(c(3,3.5))#, col=scales::alpha("black",0.3))
add_hband(c(1.31,1.3555))#, col=scales::alpha("black",0.3))
# abline(h=c(1.31), col=scales::muted("green3"))
add_band()

matplot(y=cbind(dat$GS), x=dat$YEAR, type="l", lty=1, col=c("cyan3"), ylab="Stomatal conductance\n(mol/m2/s)", xlab="Time (years)")
matlines(y=cbind(fitted(loess(dat$GS~dat$YEAR, span=60/n_year))), x=dat$YEAR, type="l", lty=1, col="black", lwd=c(1,0.5))
add_hband(c(0.16, 0.16555))#, col=scales::alpha("cyan4", 0.6))
# abline(h=c(0.16), col=scales::muted("cyan3"))
add_band()

agb = cbind(dat3$CL+dat3$CW)
matplot(y=agb, x=dat3$YEAR, type="l", lty=1, col=c("yellow4"), ylim=c(0,max(agb)), ylab="AGB\n(kgC/m2)", xlab = "Time (years)")
add_hband(c(16.9, 20.7))#, col=scales::alpha("yellow3", 0.3))
add_band()

matplot(y=cbind(dat3$CFR), x=dat3$YEAR, type="l", lty=1, col=c("brown"), ylab="C-FR\n(kgC/m2)", xlab = "Time (years)", ylim=c(0, max(dat$CFR/1e3,0.7)))
add_hband(c(0.48, 0.66))#, col=scales::alpha("brown", 0.3))
add_band()

matplot(y=cbind(dat$VCMAX), x=dat$YEAR, type="l", lty=1, col=c("green3"), ylab="Vcmax\n(umol/m2/s)", xlab="Time (years)", ylim=c(0,60))
matlines(y=cbind(fitted(loess(dat$VCMAX~dat$YEAR, span=60/n_year))), x=dat$YEAR, type="l", lty=1, col="black", lwd=c(1,0.5))
add_hband(c(20,50)) #, col=scales::muted("green4"))
add_band()

matplot(y=cbind(dat$ET, dat$PESOIL), x=dat$YEAR, type="l", lty=1, col=c("turquoise", "khaki3"), ylab="ET, PE_soil\n(mm/day)", xlab="Time (years)")
add_band()

plot_size_dist = function(){
  # matplot(y=cbind(as.numeric(colMeans(filter(dist, V1>1100 & V1<2000)[, -c(1,2)], na.rm = T)),
  #                 as.numeric(colMeans(filter(dist, V1>2100 & V1<3000)[, -c(1,2)], na.rm = T))
  #                 )*1e-2*1e4, 
  #         x=x, type="l", log="y", lty=1, col=c("black", "brown"), 
  #         xlim=c(0.01, 2), ylim=c(1e-4, 200), ylab="Density (stems/cm/ha)", xlab="Diameter (m)")
  matplot(y=cbind(as.numeric(dist_amb[gtools::mixedsort(names(dist_amb))][-1]),
                  as.numeric(dist_ele[gtools::mixedsort(names(dist_ele))][-1])
  )*1e-2*1e4, # Convert stems m-1 m-2 --> stems cm-1 ha-1
  x=x, type="l", log="y", lty=1, col=c("black", "yellow3"), 
  xlim=c(0.01, 1.2), ylim=c(1e-4, 1000), ylab="Density\n(stems/cm/ha)", xlab="Diameter (m)", las=0)
  
  abline(v=1, col=scales::alpha("red", 0.2))
  
  xobs = c(15,25,35,45,55,65,75,85,95,105)/100
  # Data for Manaus from https://link.springer.com/article/10.1007/s00442-004-1598-z 
  yobs=c(350.5221340921042,
         132.41927918860426,
         62.62503296462008,
         29.61724892214378,
         15.095996574802413,
         5.702923697662178,
         2.3219542502889836,
         1.5968055466971947,
         0.7006940913385968,
         0.5597156879584093)/10
  points(yobs~xobs, pch=20, col=scales::alpha("grey30", 0.4), cex=1.7)
}
try(plot_size_dist())

traits_obs %>% select(Leaf.LMA..g.m2., Total.BasalArea_2017.cm2.) %>% drop_na %>% with(density(x =Leaf.LMA..g.m2., weights=Total.BasalArea_2017.cm2./sum(Total.BasalArea_2017.cm2.))) %>% plot(ylim=c(0,0.02), las=0, main="", xlab="LMA", col=NA, lwd=2)
traits_obs %>% select(Leaf.LMA..g.m2., Total.BasalArea_2017.cm2.) %>% drop_na %>% with(density(x =Leaf.LMA..g.m2., weights=Total.BasalArea_2017.cm2./sum(Total.BasalArea_2017.cm2.))) %>% polygon(col=scales::alpha("grey30", 0.2), border=scales::alpha("grey30",0.4))
try(
  dat2 %>% select(YEAR, PID, BA) %>% 
    left_join(traits, by = c("PID"="SPP", "YEAR"="YEAR")) %>% 
    filter(YEAR == min(2000, max(dat2$YEAR)-1)) %>% 
    with(density(x =LMA*1000, weights=BA/sum(BA))) %>% points(col="black", type="l", lwd=1.5)
)
try(
  dat2 %>% select(YEAR, PID, BA) %>%
    left_join(traits, by = c("PID"="SPP", "YEAR"="YEAR")) %>%
    filter(YEAR == min(3000, max(dat2$YEAR)-1)) %>%
    with(density(x =LMA*1000, weights=BA/sum(BA))) %>% points(col="yellow3", type="l", lwd=1.5)
)

traits_obs %>% select(meanWoodDensity..g.cm3., Total.BasalArea_2017.cm2.) %>% drop_na %>% with(density(x =meanWoodDensity..g.cm3.*1000, weights=Total.BasalArea_2017.cm2./sum(Total.BasalArea_2017.cm2.))) %>% plot(ylim=c(0,0.005), las=0, main="", xlab="Wood density", col=NA, lwd=2)
traits_obs %>% select(meanWoodDensity..g.cm3., Total.BasalArea_2017.cm2.) %>% drop_na %>% with(density(x =meanWoodDensity..g.cm3.*1000, weights=Total.BasalArea_2017.cm2./sum(Total.BasalArea_2017.cm2.))) %>% polygon(col=scales::alpha("grey30", 0.2), border=scales::alpha("grey30",0.4))
try(
  dat2 %>% select(YEAR, PID, BA) %>% 
    left_join(traits, by = c("PID"="SPP", "YEAR"="YEAR")) %>% 
    filter(YEAR == min(2000, max(dat2$YEAR)-1)) %>% 
    with(density(x =WD, weights=BA/sum(BA))) %>% points(col="black", type="l", lwd=1.5)
)
try(
  dat2 %>% select(YEAR, PID, BA) %>%
    left_join(traits, by = c("PID"="SPP", "YEAR"="YEAR")) %>%
    filter(YEAR == min(3000, max(dat2$YEAR)-1)) %>%
    with(density(x =WD, weights=BA/sum(BA))) %>% points(col="yellow3", type="l", lwd=1.5)
)

traits_obs %>% select(Height_Max.m., Total.BasalArea_2017.cm2.) %>% drop_na %>% with(density(x =Height_Max.m., weights=Total.BasalArea_2017.cm2./sum(Total.BasalArea_2017.cm2.))) %>% plot(ylim=c(0,0.18), las=0, main="", xlab="Max. height", col=NA, lwd=2)
traits_obs %>% select(Height_Max.m., Total.BasalArea_2017.cm2.) %>% drop_na %>% with(density(x =Height_Max.m., weights=Total.BasalArea_2017.cm2./sum(Total.BasalArea_2017.cm2.))) %>% polygon(col=scales::alpha("grey30", 0.2), border=scales::alpha("grey30",0.4))
try(
  dat2 %>% select(YEAR, PID, BA) %>% 
    left_join(traits, by = c("PID"="SPP", "YEAR"="YEAR")) %>% 
    filter(YEAR == min(2000, max(dat2$YEAR)-1)) %>% 
    with(density(x =HMAT, weights=BA/sum(BA))) %>% points(col="black", type="l", lwd=1.5)
)
try(
  dat2 %>% select(YEAR, PID, BA) %>%
    left_join(traits, by = c("PID"="SPP", "YEAR"="YEAR")) %>% 
    filter(YEAR == min(3000, max(dat2$YEAR)-1)) %>% 
    with(density(x =HMAT, weights=BA/sum(BA))) %>% points(col="yellow3", type="l", lwd=1.5)
)

# traits_obs %>% select(P50..Mpa., Total.BasalArea_2017.cm2.) %>% drop_na %>% with(density(x =P50..Mpa., weights=Total.BasalArea_2017.cm2./sum(Total.BasalArea_2017.cm2.))) %>% plot(ylim=c(0,0.5), las=0, main="", xlab="Max. height", col=NA, lwd=2)
# traits_obs %>% select(P50..Mpa., Total.BasalArea_2017.cm2.) %>% drop_na %>% with(density(x =P50..Mpa., weights=Total.BasalArea_2017.cm2./sum(Total.BasalArea_2017.cm2.))) %>% polygon(col=scales::alpha("grey30", 0.2), border=scales::alpha("grey30",0.4))
# dat2 %>% select(YEAR, PID, BA) %>% 
#   filter(YEAR == min(2000, max(dat2$YEAR)-1)) %>% 
#   left_join(traits_obs %>% select(-BA), by = c("PID"="Species")) %>% 
#   drop_na %>% 
#   with(density(x =P50..Mpa., weights=BA/sum(BA))) %>% points(col="black", type="l", lwd=1.5)
# try(
#   dat2 %>% select(YEAR, PID, BA) %>% 
#     filter(YEAR == 3000) %>% 
#     left_join(traits_obs %>% select(-BA), by = c("PID"="Species")) %>% 
#     drop_na %>% 
#     with(density(x =P50..Mpa., weights=BA/sum(BA))) %>% points(col="yellow3", type="l", lwd=1.5)
# )


# 
# traits_obs %>% select(meanWoodDensity..g.cm3., Total.BasalArea_2017.cm2.) %>% drop_na %>% with(density(x =meanWoodDensity..g.cm3.*1000, weights=Total.BasalArea_2017.cm2./sum(Total.BasalArea_2017.cm2.))) %>% plot(ylim=c(0,0.005), las=0, main="", xlab="Wood density", col=NA, lwd=2)
# traits_obs %>% select(meanWoodDensity..g.cm3., Total.BasalArea_2017.cm2.) %>% drop_na %>% with(density(x =meanWoodDensity..g.cm3.*1000, weights=Total.BasalArea_2017.cm2./sum(Total.BasalArea_2017.cm2.))) %>% polygon(col=scales::alpha("grey30", 0.2), border=scales::alpha("grey30",0.4))
# dat2 %>% select(YEAR, PID, BA) %>% 
#   filter(YEAR == min(2000, max(dat2$YEAR)-1)) %>% 
#   left_join(traits_obs %>% select(-BA), by = c("PID"="Species")) %>% 
#   left_join(traits, by = c("PID"="SPP", "YEAR"="YEAR")) %>% 
#   # drop_na %>% 
#   with(density(x =meanWoodDensity..g.cm3.*1000, weights=BA/sum(BA))) %>% points(col="black", type="l", lwd=1.5)
# try(
#   dat2 %>% select(YEAR, PID, BA) %>% 
#     filter(YEAR == min(2000, max(dat2$YEAR)-1)) %>% 
#     left_join(traits_obs %>% select(-BA), by = c("PID"="Species")) %>% 
#     left_join(traits, by = c("PID"="SPP", "YEAR"="YEAR")) %>% 
#     #drop_na %>% 
#     with(density(x =WD, weights=BA/sum(BA))) %>% points(col="yellow3", type="l", lwd=1.5)
# )
# 
# 

cwm_wd = traits_obs %>% select(meanWoodDensity..g.cm3., Total.BasalArea_2017.cm2.) %>% drop_na %>% summarise(cwm=sum(meanWoodDensity..g.cm3.*Total.BasalArea_2017.cm2.)/sum(Total.BasalArea_2017.cm2.))*1000
cwm_wd_pred = dat2 %>% select(YEAR, PID, BA) %>% 
  left_join(traits, by = c("PID"="SPP", "YEAR"="YEAR")) %>% 
  drop_na %>% group_by(YEAR) %>% 
  summarize(cwm_wd = sum(WD*BA)/sum(BA))

# cwm_wd_pred %>% with(plot(cwm_wd~YEAR, type="l")) #, ylim=c(200,900)))
# add_hband(c(cwm_wd,cwm_wd+5))

# traits %>% select(YEAR, SPP, HMAT) %>% pivot_wider(names_from = "SPP", values_from = "HMAT") %>% with(matplot(x=.[,1], y=.[,2:5], lty=1, type="l"))
# traits %>% select(YEAR, SPP, WD) %>% pivot_wider(names_from = "SPP", values_from = "WD") %>% with(matplot(x=.[,1], y=.[,2:5], lty=1, type="l"))

hmat = traits %>% filter(!grepl("probe", .$SPP)) %>% select(YEAR, SPP, HMAT) %>% pivot_wider(names_from = "SPP", values_from = "HMAT") 
wd = traits %>% filter(!grepl("probe", .$SPP)) %>% select(YEAR, SPP, WD) %>% pivot_wider(names_from = "SPP", values_from = "WD") 
matplot(y=hmat[,-1], x=wd[,-1], lty=1, type="o", pch=20, cex=0.5, col=col_species, xlab="Wood density", ylab="Max height")

matplot(x=wd[,1], y=wd[,-1], col=col_species, lty=1, type="l", ylab="Wood density", xlab="Year")
matplot(x=hmat[,1], y=hmat[,-1], col=col_species, lty=1, type="l", ylab="Max. height", xlab="Year")

zz = traits %>% 
  select(YEAR, SPP, r0_avg) %>% 
  filter(!grepl(x = SPP, "probe")) %>% 
  pivot_wider(names_from = SPP, values_from = r0_avg) %>%
  as.matrix()  
zz_smooth = zz %>% as.data.frame() %>% pivot_longer(-YEAR) %>% drop_na() %>% group_by(name) %>% mutate(value = loess(value~YEAR, span=60/length(value)) %>% fitted()) %>% pivot_wider(names_from=name) %>% as.matrix

# matplot(y=tanh(zz[,-1]*20)/20, x=zz[,1], type="l", lty=1, ylab="r0", col=col_species)
matplot(y=tanh(zz_smooth[,-1]*20)/20, x=zz[,1], type="l", lty=1, ylab="r0", col=(col_species))
abline(h=0, col="black", lwd=0.2)
# abline(v=1000, col="grey")

# cwm_p50 = traits_obs %>% select(P50..Mpa., Total.BasalArea_2017.cm2.) %>% drop_na %>% summarise(cwm=sum(P50..Mpa.*Total.BasalArea_2017.cm2.)/sum(Total.BasalArea_2017.cm2.))
# cwm_p50_pred = dat2 %>% select(YEAR, PID, BA) %>% 
#   left_join(traits_obs, by = c("PID"="Species")) %>% 
#   drop_na %>% group_by(YEAR) %>% 
#   summarize(cwm_p50 = sum(P50..Mpa.*BA)/sum(BA))
# 
# cwm_p50_pred %>% with(plot(cwm_p50~YEAR, type="l")) #, ylim=c(200,900)))
# add_hband(c(cwm_wd,cwm_wd+5))


# p50 = traits %>% filter(!grepl("probe", .$SPP)) %>% select(YEAR, SPP, P50X) %>% pivot_wider(names_from = "SPP", values_from = "P50X") 
# matplot(x=p50[,1], y=p50[,-1], col=col_species, lty=1, type="l", ylab="P50", xlab="Year")
# 
# p50 %>% tail() %>% print()

if (plot_to_file) dev.off()


#### Fittest traits

#### TRAIT SPACE ####
if (plot_trait_space){

  
p1 = traits %>% 
  filter(YEAR > 50) %>% 
  # filter(RES==T) %>% 
  ggplot(aes(y=LMA, x=HMAT))+
  theme_classic(base_size = 12)+
  geom_point(aes(col=YEAR, size=RES), alpha=0.4)+
  scale_color_viridis_c(direction = -1)+
  scale_size("size_RES", range = c(0, 1.5))
# print(p1)

p1 = traits %>% 
  filter(YEAR > 50) %>% 
  filter(RES==1) %>%
  ggplot(aes(y=P50X, x=YEAR, group=SPP))+
  theme_classic(base_size = 12)+
  geom_point(aes(col=YEAR))+
  scale_color_viridis_c(direction = -1)+
  scale_size("size_RES", range = c(0, 1.5))
print(p1)
  
# matplot(y=cbind(dat$ET), x=dat$YEAR, type="l", lty=1, col=c("blue"))
  

p2 = dat2 %>% select(YEAR, PID, BA) %>% 
  left_join(traits, by = c("PID"="SPP", "YEAR"="YEAR")) %>% 
  filter(YEAR > 40 & YEAR < 50) %>% 
  # filter(RES==T) %>% 
  ggplot(aes(y=LMA, x=HMAT))+
  theme_classic(base_size = 12)+
  geom_point(aes(col=BA*1e4, size=RES), alpha=0.7)+
  scale_color_viridis_c(direction = -1)+
  scale_size("size_RES", range = c(0, 1.5), guide = F)+
  labs(col="BA")+
  ggtitle("Ambient")

p3 = dat2 %>% select(YEAR, PID, BA) %>% 
  left_join(traits, by = c("PID"="SPP", "YEAR"="YEAR")) %>% 
  filter(YEAR > 90 & YEAR < 100) %>% 
  # filter(RES==T) %>% 
  ggplot(aes(y=LMA, x=HMAT))+
  theme_classic(base_size = 12)+
  geom_point(aes(col=BA*1e4, size=RES), alpha=0.7)+
  scale_color_viridis_c(direction = -1)+
  scale_size("size_RES", range = c(0, 1.5), guide = F)+
  labs(col="BA")+
  ggtitle("Elevated")

p4 = traits_obs[1:100,] %>%  
  ggplot(aes(y=Leaf.LMA..g.m2., x=Height_Max.m.))+
  theme_classic(base_size = 12)+
  geom_point(aes(col=Total.BasalArea_2017.cm2./1e4), alpha=0.7)+
  scale_color_viridis_c(direction = -1)+
  scale_size("size_RES", range = c(0, 1.5), guide = F)+
  labs(col="BA")+
  ggtitle("Observed")

if (plot_to_file) png("traitspace.png", width=1618*1.5, height = 1196*1.5, res=300)
print(
cowplot::plot_grid(p2,p3,p4, align="hv")
)
if (plot_to_file) dev.off()

}


#### Sample results  ####
plot_sample=F

if (plot_sample){
par(mfrow=c(1,3), mar=c(5,6,4,1), oma=c(1,1,2,1), cex.lab=1.3, cex.axis=1.2, mgp=c(3.2,1,0), las=1)

with(dat %>% filter(YEAR<1200), matplot(y=cbind(GPP, NPP)*1e-3*365, x=YEAR, type="l", lty=1, col=c("green4", "green3"), ylab="GPP, NPP\n(kgC m-2 yr-1)", xlab="Time (years)"))
# abline(h=c(3,3.5), col="grey")
add_hband(c(3,3.5), col=scales::alpha("black",0.2))
add_hband(c(1.31,1.4), col=scales::alpha("black", 0.4))#, col=scales::alpha("black",0.3))
# abline(h=c(1.31), col=scales::muted("green3"))
mtext(text = "CO2 Fluxes", side=3, line=1)

traits_obs %>% select(meanWoodDensity..g.cm3., Total.BasalArea_2017.cm2.) %>% drop_na %>% with(density(x =meanWoodDensity..g.cm3.*1000, weights=Total.BasalArea_2017.cm2./sum(Total.BasalArea_2017.cm2.))) %>% plot(ylim=c(0,0.005), las=0, main="", xlab="Wood density", col=NA, lwd=2)
traits_obs %>% select(meanWoodDensity..g.cm3., Total.BasalArea_2017.cm2.) %>% drop_na %>% with(density(x =meanWoodDensity..g.cm3.*1000, weights=Total.BasalArea_2017.cm2./sum(Total.BasalArea_2017.cm2.))) %>% polygon(col=scales::alpha("grey30", 0.2), border=scales::alpha("grey30",0.4))
dat2 %>% select(YEAR, PID, BA) %>%
  left_join(traits, by = c("PID"="SPP", "YEAR"="YEAR")) %>%
  filter(YEAR == 2000) %>%
  with(density(x =WD, weights=BA/sum(BA))) %>% points(col="black", type="l", lwd=1.5)
mtext(text = "Sample\ntrait distribution", side=3, line=1)


matplot(y=cbind(as.numeric(dist_amb[gtools::mixedsort(names(dist_amb))][-1])
)*1e-2*1e4, # Convert stems m-1 m-2 --> stems cm-1 ha-1
x=x, type="l", log="y", lty=1, col=c("black", "yellow3"),
xlim=c(0.01, 1.2), ylim=c(1e-4, 1000), ylab="Density\n(stems cm-1 ha-1)", xlab="Diameter (m)", las=0)

# abline(v=1, col=scales::alpha("red", 0.2))

xobs = c(15,25,35,45,55,65,75,85,95,105)/100
# Data for Manaus from https://link.springer.com/article/10.1007/s00442-004-1598-z
yobs=c(350.5221340921042,
       132.41927918860426,
       62.62503296462008,
       29.61724892214378,
       15.095996574802413,
       5.702923697662178,
       2.3219542502889836,
       1.5968055466971947,
       0.7006940913385968,
       0.5597156879584093)/10
points(yobs~xobs, pch=20, col=scales::alpha("grey30", 0.4), cex=1.7)
mtext(text = "Size distribution", side=3, line=1)
}


## Top 5 species by BA

dat2 %>% select(YEAR, PID, BA) %>% 
  mutate(BA=BA*1e4) %>% 
  left_join(traits, by = c("PID"="SPP", "YEAR"="YEAR")) %>% 
  filter(YEAR == max(YEAR)-1) %>% arrange(desc(BA)) %>% slice(1:5)

