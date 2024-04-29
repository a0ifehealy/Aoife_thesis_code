T_ba = 1.0

#### LV chamber parameters ####
# v0_lv_ba = 5.0*1.5
v0_lv_ba = 5.0
p0_lv_ba = 1.0
Emin_lv_ba = 0.1
Emax_lv_ba = 2.5*1.75
# Emax_lv_ba = 2.5
Tes_lv_ba = 0.3
Ted_lv_ba = 0.45
Eshift_lv_ba = 0.0

#### RV Chamber parameters ####
v0_rv_ba = 10.0
p0_rv_ba = 1.0
Emin_rv_ba = 0.1
Emax_rv_ba = 1.15*1.1
# Emax_rv_ba = 1.15
Tes_rv_ba = 0.3
Ted_rv_ba = 0.45
Eshift_rv_ba = 0.0

### LA Atrium Parameters ####
# v0_la_ba = 4.0*0.75
v0_la_ba = 4.0
p0_la_ba = 1.0
Emin_la_ba = 0.15
Emax_la_ba = 0.25
Tpwb_la_ba = 0.92
Tpww_la_ba = 0.09
Tes_la_ba = Tpww_la/2
Ted_la_ba = Tpww_la
Eshift_la_ba = Tpwb_la

### RA Atrium parameters ####
v0_ra_ba = 4.0*0.75
# v0_ra_ba = 4.0
p0_ra_ba = 1.0
Emin_ra_ba = 0.15
Emax_ra_ba = 0.25
Tpwb_ra_ba = 0.92
Tpww_ra_ba = 0.09
Tes_ra_ba = Tpww_ra/2
Ted_ra_ba = Tpww_ra
Eshift_ra_ba = Tpwb_ra

#### Valve parameters ####
CQ_AV_ba = 350.0
CQ_MV_ba = 400.0
CQ_TV_ba = 400.0
CQ_PV_ba = 350.0

## Systemic Aortic Sinus ####
Csas_ba = 0.08*0.375
Rsas_ba = 0.003*16
Lsas_ba = (6.2e-5)*4
# Csas_ba = 0.08
# Rsas_ba = 0.003
# Lsas_ba = (6.2e-5)
pt0sas_ba = 100.0
qt0sas_ba = 0.0

## Systemic Artery ####
Csat_ba = 1.6
Rsat_ba = 0.05
Lsat_ba = 0.0017
pt0sat_ba = 100.0
qt0sat_ba = 0.0

## Systemic Arteriole ####
Rsar_ba = 0.5

## Systemic Capillary #### 
Rscp_ba = 0.52

## Systemic Vein ####
Csvn_ba = 20.5
Rsvn_ba = 0.075
pt0svn_ba = 0.0
qt0svn_ba = 0.0

## Pulmonary Aortic Sinus ####
Cpas_ba = 0.18*0.611
Rpas_ba = 0.002*8.352
Lpas_ba = (5.2e-5)*2.89
# Cpas_ba = 0.18
# Rpas_ba = 0.00
# Lpas_ba = (5.2e-5)
pt0pas_ba = 30.0
qt0pas_ba = 0.0

## Pulmonary Artery ####
Cpat_ba = 3.8
Rpat_ba = 0.01
Lpat_ba = 0.0017
pt0pat_ba = 30.0
qt0pat_ba = 0.0

## Pulmonary Arteriole ####
Rpar_ba = 0.05

## Pulmonary Capillary ####
Rpcp_ba = 0.25

## Pulmonary Vein ####
Cpvn_ba = 20.5
Rpvn_ba = 0.006 
pt0pvn_ba = 0.0
qt0pvn_ba = 0.0

# Valve initial conditions #### 
LV_Vt0_ba = 500
RV_Vt0_ba = 500
LA_Vt0_ba = 20
RA_Vt0_ba = 20