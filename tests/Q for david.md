Q for david.

0) bug with AI
1) sand+clay+silt add up to 1, so what about OM and gravel? 
2) psi_m should be re-calculated at the end of the step?
3) is n starting with 0 or 1 - I get correct results with 1
4) shouldnt solar calcs be done at the start of the step rather than in the middle of the step - latter implies that half of input vars are updated and half arent

Changes:
1) saved end of timestep stress_factor = (sm-RES)/(Wmax-RES) and psi_m = f(theta_i, etc)
2) used this stress_factor to do solar Rad calcs in beginning of next timestep
3) saved ppfd to srad struct (for use in T calcs)
4) modified splash.calc_daily_fluxes to return etr, so that aet can be returned without breaking current code
4) modified evap.calc_daily_fluxes to input T, and fed it sw_in*(1-fapar) instead of sw_in
