T = 1.0

#### Windkessel Parameters ####
wk_rp = 0.9
wk_rc = 0.05
wk_c = 1.1
wk_l = 0.0025
wk_A_sine = 90
wk_A_trunc = 90*pi
wk_tspan = (0, 10.0)
wk_capacitor_delp0 = -78
wk_reindex_t = 7

#### Systemic Pulmonary Parameters ####
A_sys_pul = 90*pi
sys_pul_tspan = (0, 20.0)
sys_pul_reindex_t = 17

#### LV chamber parameters ####
v0_lv = 5.0
p0_lv = 1.0
Emin_lv = 0.1
Emax_lv = 2.5
Tes_lv = 0.3
Ted_lv = 0.45
Eshift_lv = 0.0

#### RV Chamber parameters ####
v0_rv = 10.0
p0_rv = 1.0
Emin_rv = 0.1
Emax_rv = 1.15
Tes_rv = 0.3
Ted_rv = 0.45
Eshift_rv = 0.0

### LA Atrium Parameters ####
v0_la = 4.0
p0_la = 1.0
Emin_la = 0.15
Emax_la = 0.25
Tpwb_la = 0.92
Tpww_la = 0.09
Tes_la = Tpww_la/2
Ted_la = Tpww_la
Eshift_la = Tpwb_la

### RA Atrium parameters ####
v0_ra = 4.0
p0_ra = 1.0
Emin_ra = 0.15
Emax_ra = 0.25
Tpwb_ra = 0.92
Tpww_ra = 0.09
Tes_ra = Tpww_ra/2
Ted_ra = Tpww_ra
Eshift_ra = Tpwb_ra

#### Valve parameters ####
CQ_AV = 350.0
CQ_MV = 400.0
CQ_TV = 400.0
CQ_PV = 350.0

## Systemic Aortic Sinus ####
Csas = 0.08
Rsas = 0.003
Lsas = 6.2e-5
pt0sas = 100.0
qt0sas = 0.0

## Systemic Artery ####
Csat = 1.6
Rsat = 0.05
Lsat = 0.0017
pt0sat = 100.0
qt0sat = 0.0

## Systemic Arteriole ####
Rsar = 0.5

## Systemic Capillary #### 
Rscp = 0.52

## Systemic Vein ####
Csvn = 20.5
Rsvn = 0.075
pt0svn = 0.0
qt0svn = 0.0

## Pulmonary Aortic Sinus ####
Cpas = 0.18
Rpas = 0.002
Lpas = 5.2e-5
pt0pas = 30.0
qt0pas = 0.0

## Pulmonary Artery ####
Cpat = 3.8
Rpat = 0.01
Lpat = 0.0017
pt0pat = 30.0
qt0pat = 0.0

## Pulmonary Arteriole ####
Rpar = 0.05

## Pulmonary Capillary ####
Rpcp = 0.25

## Pulmonary Vein ####
Cpvn = 20.5
Rpvn = 0.006 
pt0pvn = 0.0
qt0pvn = 0.0

# Valve initial conditions #### 
LV_Vt0 = 500
RV_Vt0 = 500
LA_Vt0 = 20
RA_Vt0 = 20