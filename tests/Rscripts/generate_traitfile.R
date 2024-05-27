library(dplyr)

n = 5
filename = "~/codes/Plant-FATE/tests/data/traits_droughtmip.csv"

# Default traits
tibble(species = paste0("spp", 1:n)) |>
  mutate (
      lma          = 0.119        # leaf mass per leaf area [kg/m2]
    , zeta         = 0.2		    	# fine root mass per leaf area [kg/m2] 0.7
    , fcr          = 0.47         # coarse root mass per unit stem mass   [REF: GRoot database]
    , hmat         = 23.99		    # maximum height [m]
    , fhmat        = 0.8          # height at reproductive maturity as fraction of hmat [Guesstimate]
    , seed_mass    = 3.8e-5	      # [kg]
    , wood_density = 690	        # [kg/m3]
    , p50_xylem    = -2.29        # Xylem P50 [MPa]  -4.515
    , K_leaf       = 0.5e-16		  # Leaf conductance [m]  ---> ** Calibrated to gs ** 
    , K_xylem      = 4e-16		    # Leaf conductance [m]
    , b_leaf       = 1			      # Shape parameter of leaf vulnerabilty curve [-] (Ref: Joshi et al 2022)
    , b_xylem      = 1            # Shape parameter of xylem vulnerabilty curve [-] 
    , sm_xylem     = 0            # Xylem safety margin (P50x - Pg88) [MPa]
    , a            = 75           # height-diameter allometry 114
    , c            = 6000         # crown area allometry
    , m            = 1.5	        # crown shape smoothness
    , n            = 3            # crown top-heaviness
    # , fg           = 0.1          # upper canopy gap fraction
  ) |> 
  # mutate(lma = runif(n, min=40, max=270)) |>
  mutate(hmat = runif(n, min=10, max=30)) |>
  mutate(wood_density = runif(n, min=300, max=1200)) |>
  mutate(p50_xylem = runif(n, min=-5, max=-0.5)) |>
  mutate(sm_xylem = runif(n, min=-1, max=-p50_xylem)) |>
  write.csv(file=filename, row.names = F)
  

