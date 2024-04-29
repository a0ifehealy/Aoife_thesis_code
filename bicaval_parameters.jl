T_bc = 1.0

#### LV chamber parameters ####
v0_lv_bc = 5.0*1.5
p0_lv_bc = 1.0
Emin_lv_bc = 0.1
Emax_lv_bc = 2.5*1.75
Tes_lv_bc = 0.3
Ted_lv_bc = 0.45
Eshift_lv_bc = 0.0

#### RV Chamber parameters ####
v0_rv_bc = 10.0
p0_rv_bc = 1.0
Emin_rv_bc = 0.1
Emax_rv_bc = 1.15*1.1
Tes_rv_bc = 0.3
Ted_rv_bc = 0.45
Eshift_rv_bc = 0.0

### LA Atrium Parameters ####
v0_la_bc = 4.0*0.75
p0_la_bc = 1.0
Emin_la_bc = 0.15
Emax_la_bc = 0.25
Tpwb_la_bc = 0.92
Tpww_la_bc = 0.09
Tes_la_bc = Tpww_la/2
Ted_la_bc = Tpww_la
Eshift_la_bc = Tpwb_la

### RA Atrium parameters ####
v0_ra_bc = 4.0*0.5
p0_ra_bc = 1.0
Emin_ra_bc = 0.15
Emax_ra_bc = 0.25
Tpwb_ra_bc = 0.92
Tpww_ra_bc = 0.09
Tes_ra_bc = Tpww_ra/2
Ted_ra_bc = Tpww_ra
Eshift_ra_bc = Tpwb_ra

#### Valve parameters ####
CQ_AV_bc = 350.0
CQ_MV_bc = 400.0
CQ_TV_bc = 400.0
CQ_PV_bc = 350.0

## Systemic Aortic Sinus ####
Csas_bc = 0.08*0.375
Rsas_bc = 0.003*16
Lsas_bc = (6.2e-5)*4
pt0sas_bc = 100.0
qt0sas_bc = 0.0

## Systemic Artery ####
Csat_bc = 1.6
Rsat_bc = 0.05
Lsat_bc = 0.0017
pt0sat_bc = 100.0
qt0sat_bc = 0.0

## Systemic Arteriole ####
Rsar_bc = 0.5

## Systemic Capillary #### 
Rscp_bc = 0.52

## Systemic Vein ####
Csvn_bc = (20.5*0.95)+(20.5*0.05*3*0.614)
Rsvn_bc = (0.075*0.995)+(0.075*0.005*1.916)
pt0svn_bc = 0.0
qt0svn_bc = 0.0

## Pulmonary Aortic Sinus ####
Cpas_bc = 0.18*0.611
Rpas_bc = 0.002*8.352
Lpas_bc = (5.2e-5)*2.89
pt0pas_bc = 30.0
qt0pas_bc = 0.0

## Pulmonary Artery ####
Cpat_bc = 3.8
Rpat_bc = 0.01
Lpat_bc = 0.0017
pt0pat_bc = 30.0
qt0pat_bc = 0.0

## Pulmonary Arteriole ####
Rpar_bc = 0.05

## Pulmonary Capillary ####
Rpcp_bc = 0.25

## Pulmonary Vein ####
Cpvn_bc = 20.5
Rpvn_bc = 0.006 
pt0pvn_bc = 0.0
qt0pvn_bc = 0.0

# Valve initial conditions #### 
LV_Vt0_bc = 500
RV_Vt0_bc = 500
LA_Vt0_bc = 20
RA_Vt0_bc = 20