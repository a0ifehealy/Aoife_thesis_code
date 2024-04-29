#############################################################################
# This script outputs graphs for different parts of the circulatory system 
# with 12 sub graphsy showing the effect of 12 different physical changes in 
# isolation
#############################################################################

using ModelingToolkit
using Plots
using DifferentialEquations

@variables t
D = Differential(t)
include("healthy_human_parameters.jl")
include("circulatory_system_components.jl")


######################
# ELASTANCE FUNCTION #
######################

function Elastance(t, E_min, E_max, T, T_es, T_ep, Eshift)
    # Get time relative to start of cycle
    t_i = rem(t + (1 - Eshift) * T, T)

    E_p = (t_i <= T_es) * (1 - cos(t_i / T_es * pi)) / 2 +
         (t_i > T_es) * (t_i <= T_ep) * (1 + cos((t_i - T_es) / (T_ep - T_es) * pi)) / 2 +
         (t_i > T_ep) * 0

    E = E_min + (E_max - E_min) * E_p

    return E
end


##########################
# INITIALISING ALL PLOTS #
##########################

function plot_setup(plot_title, x_label, y_label)
    return plot(xlabel=x_label, ylabel=y_label, legend=false, title=plot_title, yguidefontsize=12, ytickfontsize=12, xguidefontsize=12, xtickfontsize=12)
end

plot_title = "(a) Decreased aortic sinus\nradius"

p_rv_pt_ao_r = plot_setup(plot_title, "t (s)", "Pressure (mmHg)")
p_ra_pt_ao_r = plot_setup(plot_title, "t (s)", "Pressure (mmHg)")
p_lv_pt_ao_r = plot_setup(plot_title, "t (s)", "Pressure (mmHg)")
p_la_pt_ao_r = plot_setup(plot_title, "t (s)", "Pressure (mmHg)")
p_rv_pv_ao_r = plot_setup(plot_title, "Volume (ml)", "Pressure (mmHg)")
p_ra_pv_ao_r = plot_setup(plot_title, "Volume (ml)", "Pressure (mmHg)")
p_lv_pv_ao_r = plot_setup(plot_title, "Volume (ml)", "Pressure (mmHg)")
p_la_pv_ao_r = plot_setup(plot_title, "Volume (ml)", "Pressure (mmHg)")
p_ao_pt_ao_r = plot_setup(plot_title, "t (s)", "Pressure (mmHg)")
p_pa_pt_ao_r = plot_setup(plot_title, "t (s)", "Pressure (mmHg)")
p_lv_vt_ao_r = plot_setup(plot_title, "t (s)", "Volume (ml)")
p_rv_vt_ao_r = plot_setup(plot_title, "t (s)", "Volume (ml)")
p_la_vt_ao_r = plot_setup(plot_title, "t (s)", "Volume (ml)")
p_ra_vt_ao_r = plot_setup(plot_title, "t (s)", "Volume (ml)")
p_ao_vt_ao_r = plot_setup(plot_title, "t (s)", "Volume (ml)")
p_pa_vt_ao_r = plot_setup(plot_title, "t (s)", "Volume (ml)");

plot_title = "(b) Decreased aortic sinus\nYoung's modulus"

p_rv_pt_ao_ym = plot_setup(plot_title, "t (s)", "Pressure (mmHg)")
p_ra_pt_ao_ym = plot_setup(plot_title, "t (s)", "Pressure (mmHg)")
p_lv_pt_ao_ym = plot_setup(plot_title, "t (s)", "Pressure (mmHg)")
p_la_pt_ao_ym = plot_setup(plot_title, "t (s)", "Pressure (mmHg)")
p_rv_pv_ao_ym = plot_setup(plot_title, "Volume (ml)", "Pressure (mmHg)")
p_ra_pv_ao_ym = plot_setup(plot_title, "Volume (ml)", "Pressure (mmHg)")
p_lv_pv_ao_ym = plot_setup(plot_title, "Volume (ml)", "Pressure (mmHg)")
p_la_pv_ao_ym = plot_setup(plot_title, "Volume (ml)", "Pressure (mmHg)")
p_ao_pt_ao_ym = plot_setup(plot_title, "t (s)", "Pressure (mmHg)")
p_pa_pt_ao_ym = plot_setup(plot_title, "t (s)", "Pressure (mmHg)")
p_lv_vt_ao_ym = plot_setup(plot_title, "t (s)", "Volume (ml)")
p_rv_vt_ao_ym = plot_setup(plot_title, "t (s)", "Volume (ml)")
p_la_vt_ao_ym = plot_setup(plot_title, "t (s)", "Volume (ml)")
p_ra_vt_ao_ym = plot_setup(plot_title, "t (s)", "Volume (ml)")
p_ao_vt_ao_ym = plot_setup(plot_title, "t (s)", "Volume (ml)")
p_pa_vt_ao_ym = plot_setup(plot_title, "t (s)", "Volume (ml)");

plot_title = "(c) Decreased PA radius"

p_rv_pt_pa_r = plot_setup(plot_title, "t (s)", "Pressure (mmHg)")
p_ra_pt_pa_r = plot_setup(plot_title, "t (s)", "Pressure (mmHg)")
p_lv_pt_pa_r = plot_setup(plot_title, "t (s)", "Pressure (mmHg)")
p_la_pt_pa_r = plot_setup(plot_title, "t (s)", "Pressure (mmHg)")
p_rv_pv_pa_r = plot_setup(plot_title, "Volume (ml)", "Pressure (mmHg)")
p_ra_pv_pa_r = plot_setup(plot_title, "Volume (ml)", "Pressure (mmHg)")
p_lv_pv_pa_r = plot_setup(plot_title, "Volume (ml)", "Pressure (mmHg)")
p_la_pv_pa_r = plot_setup(plot_title, "Volume (ml)", "Pressure (mmHg)")
p_ao_pt_pa_r = plot_setup(plot_title, "t (s)", "Pressure (mmHg)")
p_pa_pt_pa_r = plot_setup(plot_title, "t (s)", "Pressure (mmHg)")
p_lv_vt_pa_r = plot_setup(plot_title, "t (s)", "Volume (ml)")
p_rv_vt_pa_r = plot_setup(plot_title, "t (s)", "Volume (ml)")
p_la_vt_pa_r = plot_setup(plot_title, "t (s)", "Volume (ml)")
p_ra_vt_pa_r = plot_setup(plot_title, "t (s)", "Volume (ml)")
p_ao_vt_pa_r = plot_setup(plot_title, "t (s)", "Volume (ml)")
p_pa_vt_pa_r = plot_setup(plot_title, "t (s)", "Volume (ml)");

plot_title = "(d) Decreased PA Young's\nmodulus"

p_rv_pt_pa_ym = plot_setup(plot_title, "t (s)", "Pressure (mmHg)")
p_ra_pt_pa_ym = plot_setup(plot_title, "t (s)", "Pressure (mmHg)")
p_lv_pt_pa_ym = plot_setup(plot_title, "t (s)", "Pressure (mmHg)")
p_la_pt_pa_ym = plot_setup(plot_title, "t (s)", "Pressure (mmHg)")
p_rv_pv_pa_ym = plot_setup(plot_title, "Volume (ml)", "Pressure (mmHg)")
p_ra_pv_pa_ym = plot_setup(plot_title, "Volume (ml)", "Pressure (mmHg)")
p_lv_pv_pa_ym = plot_setup(plot_title, "Volume (ml)", "Pressure (mmHg)")
p_la_pv_pa_ym = plot_setup(plot_title, "Volume (ml)", "Pressure (mmHg)")
p_ao_pt_pa_ym = plot_setup(plot_title, "t (s)", "Pressure (mmHg)")
p_pa_pt_pa_ym = plot_setup(plot_title, "t (s)", "Pressure (mmHg)")
p_lv_vt_pa_ym = plot_setup(plot_title, "t (s)", "Volume (ml)")
p_rv_vt_pa_ym = plot_setup(plot_title, "t (s)", "Volume (ml)")
p_la_vt_pa_ym = plot_setup(plot_title, "t (s)", "Volume (ml)")
p_ra_vt_pa_ym = plot_setup(plot_title, "t (s)", "Volume (ml)")
p_ao_vt_pa_ym = plot_setup(plot_title, "t (s)", "Volume (ml)")
p_pa_vt_pa_ym = plot_setup(plot_title, "t (s)", "Volume (ml)");

plot_title = "(e) Decreased LA volume"

p_rv_pt_la_v0 = plot_setup(plot_title, "t (s)", "Pressure (mmHg)")
p_ra_pt_la_v0 = plot_setup(plot_title, "t (s)", "Pressure (mmHg)")
p_lv_pt_la_v0 = plot_setup(plot_title, "t (s)", "Pressure (mmHg)")
p_la_pt_la_v0 = plot_setup(plot_title, "t (s)", "Pressure (mmHg)")
p_rv_pv_la_v0 = plot_setup(plot_title, "Volume (ml)", "Pressure (mmHg)")
p_ra_pv_la_v0 = plot_setup(plot_title, "Volume (ml)", "Pressure (mmHg)")
p_lv_pv_la_v0 = plot_setup(plot_title, "Volume (ml)", "Pressure (mmHg)")
p_la_pv_la_v0 = plot_setup(plot_title, "Volume (ml)", "Pressure (mmHg)")
p_ao_pt_la_v0 = plot_setup(plot_title, "t (s)", "Pressure (mmHg)")
p_pa_pt_la_v0 = plot_setup(plot_title, "t (s)", "Pressure (mmHg)")
p_lv_vt_la_v0 = plot_setup(plot_title, "t (s)", "Volume (ml)")
p_rv_vt_la_v0 = plot_setup(plot_title, "t (s)", "Volume (ml)")
p_la_vt_la_v0 = plot_setup(plot_title, "t (s)", "Volume (ml)")
p_ra_vt_la_v0 = plot_setup(plot_title, "t (s)", "Volume (ml)")
p_ao_vt_la_v0 = plot_setup(plot_title, "t (s)", "Volume (ml)")
p_pa_vt_la_v0 = plot_setup(plot_title, "t (s)", "Volume (ml)");

plot_title = "(f) Decreased RA volume\n(biatrial)"

p_rv_pt_ra_v0_ba = plot_setup(plot_title, "t (s)", "Pressure (mmHg)")
p_ra_pt_ra_v0_ba = plot_setup(plot_title, "t (s)", "Pressure (mmHg)")
p_lv_pt_ra_v0_ba = plot_setup(plot_title, "t (s)", "Pressure (mmHg)")
p_la_pt_ra_v0_ba = plot_setup(plot_title, "t (s)", "Pressure (mmHg)")
p_rv_pv_ra_v0_ba = plot_setup(plot_title, "Volume (ml)", "Pressure (mmHg)")
p_ra_pv_ra_v0_ba = plot_setup(plot_title, "Volume (ml)", "Pressure (mmHg)")
p_lv_pv_ra_v0_ba = plot_setup(plot_title, "Volume (ml)", "Pressure (mmHg)")
p_la_pv_ra_v0_ba = plot_setup(plot_title, "Volume (ml)", "Pressure (mmHg)")
p_ao_pt_ra_v0_ba = plot_setup(plot_title, "t (s)", "Pressure (mmHg)")
p_pa_pt_ra_v0_ba = plot_setup(plot_title, "t (s)", "Pressure (mmHg)")
p_lv_vt_ra_v0_ba = plot_setup(plot_title, "t (s)", "Volume (ml)")
p_rv_vt_ra_v0_ba = plot_setup(plot_title, "t (s)", "Volume (ml)")
p_la_vt_ra_v0_ba = plot_setup(plot_title, "t (s)", "Volume (ml)")
p_ra_vt_ra_v0_ba = plot_setup(plot_title, "t (s)", "Volume (ml)")
p_ao_vt_ra_v0_ba = plot_setup(plot_title, "t (s)", "Volume (ml)")
p_pa_vt_ra_v0_ba = plot_setup(plot_title, "t (s)", "Volume (ml)");

plot_title = "(g) Increased LV elastance"

p_rv_pt_lv_emax = plot_setup(plot_title, "t (s)", "Pressure (mmHg)")
p_ra_pt_lv_emax = plot_setup(plot_title, "t (s)", "Pressure (mmHg)")
p_lv_pt_lv_emax = plot_setup(plot_title, "t (s)", "Pressure (mmHg)")
p_la_pt_lv_emax = plot_setup(plot_title, "t (s)", "Pressure (mmHg)")
p_rv_pv_lv_emax = plot_setup(plot_title, "Volume (ml)", "Pressure (mmHg)")
p_ra_pv_lv_emax = plot_setup(plot_title, "Volume (ml)", "Pressure (mmHg)")
p_lv_pv_lv_emax = plot_setup(plot_title, "Volume (ml)", "Pressure (mmHg)")
p_la_pv_lv_emax = plot_setup(plot_title, "Volume (ml)", "Pressure (mmHg)")
p_ao_pt_lv_emax = plot_setup(plot_title, "t (s)", "Pressure (mmHg)")
p_pa_pt_lv_emax = plot_setup(plot_title, "t (s)", "Pressure (mmHg)")
p_lv_vt_lv_emax = plot_setup(plot_title, "t (s)", "Volume (ml)")
p_rv_vt_lv_emax = plot_setup(plot_title, "t (s)", "Volume (ml)")
p_la_vt_lv_emax = plot_setup(plot_title, "t (s)", "Volume (ml)")
p_ra_vt_lv_emax = plot_setup(plot_title, "t (s)", "Volume (ml)")
p_ao_vt_lv_emax = plot_setup(plot_title, "t (s)", "Volume (ml)")
p_pa_vt_lv_emax = plot_setup(plot_title, "t (s)", "Volume (ml)");

plot_title = "(h) Increased RV elastance"

p_rv_pt_rv_emax = plot_setup(plot_title, "t (s)", "Pressure (mmHg)")
p_ra_pt_rv_emax = plot_setup(plot_title, "t (s)", "Pressure (mmHg)")
p_lv_pt_rv_emax = plot_setup(plot_title, "t (s)", "Pressure (mmHg)")
p_la_pt_rv_emax = plot_setup(plot_title, "t (s)", "Pressure (mmHg)")
p_rv_pv_rv_emax = plot_setup(plot_title, "Volume (ml)", "Pressure (mmHg)")
p_ra_pv_rv_emax = plot_setup(plot_title, "Volume (ml)", "Pressure (mmHg)")
p_lv_pv_rv_emax = plot_setup(plot_title, "Volume (ml)", "Pressure (mmHg)")
p_la_pv_rv_emax = plot_setup(plot_title, "Volume (ml)", "Pressure (mmHg)")
p_ao_pt_rv_emax = plot_setup(plot_title, "t (s)", "Pressure (mmHg)")
p_pa_pt_rv_emax = plot_setup(plot_title, "t (s)", "Pressure (mmHg)")
p_lv_vt_rv_emax = plot_setup(plot_title, "t (s)", "Volume (ml)")
p_rv_vt_rv_emax = plot_setup(plot_title, "t (s)", "Volume (ml)")
p_la_vt_rv_emax = plot_setup(plot_title, "t (s)", "Volume (ml)")
p_ra_vt_rv_emax = plot_setup(plot_title, "t (s)", "Volume (ml)")
p_ao_vt_rv_emax = plot_setup(plot_title, "t (s)", "Volume (ml)")
p_pa_vt_rv_emax = plot_setup(plot_title, "t (s)", "Volume (ml)");

plot_title = "(i) Increased LV volume"

p_rv_pt_lv_v0 = plot_setup(plot_title, "t (s)", "Pressure (mmHg)")
p_ra_pt_lv_v0 = plot_setup(plot_title, "t (s)", "Pressure (mmHg)")
p_lv_pt_lv_v0 = plot_setup(plot_title, "t (s)", "Pressure (mmHg)")
p_la_pt_lv_v0 = plot_setup(plot_title, "t (s)", "Pressure (mmHg)")
p_rv_pv_lv_v0 = plot_setup(plot_title, "Volume (ml)", "Pressure (mmHg)")
p_ra_pv_lv_v0 = plot_setup(plot_title, "Volume (ml)", "Pressure (mmHg)")
p_lv_pv_lv_v0 = plot_setup(plot_title, "Volume (ml)", "Pressure (mmHg)")
p_la_pv_lv_v0 = plot_setup(plot_title, "Volume (ml)", "Pressure (mmHg)")
p_ao_pt_lv_v0 = plot_setup(plot_title, "t (s)", "Pressure (mmHg)")
p_pa_pt_lv_v0 = plot_setup(plot_title, "t (s)", "Pressure (mmHg)")
p_lv_vt_lv_v0 = plot_setup(plot_title, "t (s)", "Volume (ml)")
p_rv_vt_lv_v0 = plot_setup(plot_title, "t (s)", "Volume (ml)")
p_la_vt_lv_v0 = plot_setup(plot_title, "t (s)", "Volume (ml)")
p_ra_vt_lv_v0 = plot_setup(plot_title, "t (s)", "Volume (ml)")
p_ao_vt_lv_v0 = plot_setup(plot_title, "t (s)", "Volume (ml)")
p_pa_vt_lv_v0 = plot_setup(plot_title, "t (s)", "Volume (ml)");

plot_title = "(j) Decreased RA volume\n(bicaval)"

p_rv_pt_ra_v0_bc = plot_setup(plot_title, "t (s)", "Pressure (mmHg)")
p_ra_pt_ra_v0_bc = plot_setup(plot_title, "t (s)", "Pressure (mmHg)")
p_lv_pt_ra_v0_bc = plot_setup(plot_title, "t (s)", "Pressure (mmHg)")
p_la_pt_ra_v0_bc = plot_setup(plot_title, "t (s)", "Pressure (mmHg)")
p_rv_pv_ra_v0_bc = plot_setup(plot_title, "Volume (ml)", "Pressure (mmHg)")
p_ra_pv_ra_v0_bc = plot_setup(plot_title, "Volume (ml)", "Pressure (mmHg)")
p_lv_pv_ra_v0_bc = plot_setup(plot_title, "Volume (ml)", "Pressure (mmHg)")
p_la_pv_ra_v0_bc = plot_setup(plot_title, "Volume (ml)", "Pressure (mmHg)")
p_ao_pt_ra_v0_bc = plot_setup(plot_title, "t (s)", "Pressure (mmHg)")
p_pa_pt_ra_v0_bc = plot_setup(plot_title, "t (s)", "Pressure (mmHg)")
p_lv_vt_ra_v0_bc = plot_setup(plot_title, "t (s)", "Volume (ml)")
p_rv_vt_ra_v0_bc = plot_setup(plot_title, "t (s)", "Volume (ml)")
p_la_vt_ra_v0_bc = plot_setup(plot_title, "t (s)", "Volume (ml)")
p_ra_vt_ra_v0_bc = plot_setup(plot_title, "t (s)", "Volume (ml)")
p_ao_vt_ra_v0_bc = plot_setup(plot_title, "t (s)", "Volume (ml)")
p_pa_vt_ra_v0_bc = plot_setup(plot_title, "t (s)", "Volume (ml)");

plot_title = "(k) Decreased VC radius"

p_rv_pt_vc_r = plot_setup(plot_title, "t (s)", "Pressure (mmHg)")
p_ra_pt_vc_r = plot_setup(plot_title, "t (s)", "Pressure (mmHg)")
p_lv_pt_vc_r = plot_setup(plot_title, "t (s)", "Pressure (mmHg)")
p_la_pt_vc_r = plot_setup(plot_title, "t (s)", "Pressure (mmHg)")
p_rv_pv_vc_r = plot_setup(plot_title, "Volume (ml)", "Pressure (mmHg)")
p_ra_pv_vc_r = plot_setup(plot_title, "Volume (ml)", "Pressure (mmHg)")
p_lv_pv_vc_r = plot_setup(plot_title, "Volume (ml)", "Pressure (mmHg)")
p_la_pv_vc_r = plot_setup(plot_title, "Volume (ml)", "Pressure (mmHg)")
p_ao_pt_vc_r = plot_setup(plot_title, "t (s)", "Pressure (mmHg)")
p_pa_pt_vc_r = plot_setup(plot_title, "t (s)", "Pressure (mmHg)")
p_lv_vt_vc_r = plot_setup(plot_title, "t (s)", "Volume (ml)")
p_rv_vt_vc_r = plot_setup(plot_title, "t (s)", "Volume (ml)")
p_la_vt_vc_r = plot_setup(plot_title, "t (s)", "Volume (ml)")
p_ra_vt_vc_r = plot_setup(plot_title, "t (s)", "Volume (ml)")
p_ao_vt_vc_r = plot_setup(plot_title, "t (s)", "Volume (ml)")
p_pa_vt_vc_r = plot_setup(plot_title, "t (s)", "Volume (ml)");

plot_title = "(l) Decreased VC Young's\nmodulus"

p_rv_pt_vc_ym = plot_setup(plot_title, "t (s)", "Pressure (mmHg)")
p_ra_pt_vc_ym = plot_setup(plot_title, "t (s)", "Pressure (mmHg)")
p_lv_pt_vc_ym = plot_setup(plot_title, "t (s)", "Pressure (mmHg)")
p_la_pt_vc_ym = plot_setup(plot_title, "t (s)", "Pressure (mmHg)")
p_rv_pv_vc_ym = plot_setup(plot_title, "Volume (ml)", "Pressure (mmHg)")
p_ra_pv_vc_ym = plot_setup(plot_title, "Volume (ml)", "Pressure (mmHg)")
p_lv_pv_vc_ym = plot_setup(plot_title, "Volume (ml)", "Pressure (mmHg)")
p_la_pv_vc_ym = plot_setup(plot_title, "Volume (ml)", "Pressure (mmHg)")
p_ao_pt_vc_ym = plot_setup(plot_title, "t (s)", "Pressure (mmHg)")
p_pa_pt_vc_ym = plot_setup(plot_title, "t (s)", "Pressure (mmHg)")
p_lv_vt_vc_ym = plot_setup(plot_title, "t (s)", "Volume (ml)")
p_rv_vt_vc_ym = plot_setup(plot_title, "t (s)", "Volume (ml)")
p_la_vt_vc_ym = plot_setup(plot_title, "t (s)", "Volume (ml)")
p_ra_vt_vc_ym = plot_setup(plot_title, "t (s)", "Volume (ml)")
p_ao_vt_vc_ym = plot_setup(plot_title, "t (s)", "Volume (ml)")
p_pa_vt_vc_ym = plot_setup(plot_title, "t (s)", "Volume (ml)");


###########################
# REFERENCE HEALTHY HUMAN #
###########################

# Heart chambers
@named LV = HeartChamber(V_0=v0_lv, p_0=p0_lv, E_min=Emin_lv, E_max=Emax_lv, T=T, T_es=Tes_lv, T_ep=Ted_lv, Eshift=0.0)
@named LA = HeartChamber(V_0=v0_la, p_0=p0_la, E_min=Emin_la, E_max=Emax_la, T=T, T_es=Tpww_la / 2, T_ep=Tpww_la, Eshift=Tpwb_la)
@named RV = HeartChamber(V_0=v0_rv, p_0=p0_rv, E_min=Emin_rv, E_max=Emax_rv, T=T, T_es=Tes_rv, T_ep=Ted_rv, Eshift=0.0)
@named RA = HeartChamber(V_0=v0_ra, p_0=p0_ra, E_min=Emin_ra, E_max=Emax_ra, T=T, T_es=Tpww_ra / 2, T_ep=Tpww_ra, Eshift=Tpwb_ra)

# Valves
@named AV = HeartValve(CQ=CQ_AV)
@named MV = HeartValve(CQ=CQ_MV)
@named TV = HeartValve(CQ=CQ_TV)
@named PV = HeartValve(CQ=CQ_PV)

# Systemic circulation
@named SAS = CRL(C=Csas, R=Rsas, L=Lsas)
@named SAT = CRL(C=Csat, R=Rsat, L=Lsat)
@named SAR = Resistor(R=Rsar)
@named SCP = Resistor(R=Rscp)
@named SVN = CR(R=Rsvn, C=Csvn)

# Pulmonary circulatioin
@named PAS = CRL(C=Cpas, R=Rpas, L=Lpas)
@named PAT = CRL(C=Cpat, R=Rpat, L=Lpat)
@named PAR = Resistor(R=Rpar)
@named PCP = Resistor(R=Rpcp)
@named PVN = CR(R=Rpvn, C=Cpvn);

# Connect
circ_eqs = [
    connect(LV.out, AV.in)
    connect(AV.out, SAS.in)
    connect(SAS.out, SAT.in)
    connect(SAT.out, SAR.in)
    connect(SAR.out, SCP.in)
    connect(SCP.out, SVN.in)
    connect(SVN.out, RA.in)
    connect(RA.out, TV.in)
    connect(TV.out, RV.in)
    connect(RV.out, PV.in)
    connect(PV.out, PAS.in)
    connect(PAS.out, PAT.in)
    connect(PAT.out, PAR.in)
    connect(PAR.out, PCP.in)
    connect(PCP.out, PVN.in)
    connect(PVN.out, LA.in)
    connect(LA.out, MV.in)
    connect(MV.out, LV.in)
];

# Create ODE system
@named _circ_model = ODESystem(circ_eqs, t)
@named circ_model = compose(_circ_model,
[LV, RV, LA, RA, AV, MV, PV, TV, SAS, SAT, SAR, SCP, SVN, PAS, PAT, PAR, PCP, PVN])
circ_sys = structural_simplify(circ_model)

#Initial conditions for ODEs
u0 = [
    LV.V => LV_Vt0
    LV.p => (LV_Vt0 - v0_lv) * Emin_lv + p0_lv
    RV.V => RV_Vt0
    RV.p => (RV_Vt0 - v0_rv) * Emin_rv + p0_rv
    LA.V => LA_Vt0
    RA.V => RA_Vt0
    SAS.C.p => pt0sas
    SAS.C.V => pt0sas * Csas
    SAS.L.q => qt0sas
    SAT.C.p => pt0sat
    SAT.C.V => pt0sat * Csat
    SAT.L.q => qt0sat
    SVN.C.p => pt0svn
    SVN.C.V => pt0svn * Csvn
    PAS.C.p => pt0pas
    PAS.C.V => pt0pas * Cpas
    PAS.L.q => qt0pas
    PAT.C.p => pt0pat
    PAT.C.V => pt0pat * Cpat
    PAT.L.q => qt0pat
    PVN.C.p => pt0pvn
    PVN.C.V => pt0pvn * Cpvn
];

# Solve system of ODEs
prob = ODAEProblem(circ_sys, u0, (0.0, 20.0))
sol = solve(prob, RK4(), reltol=1e-9, abstol=1e-12, saveat=19:0.01:20);


#########################
# ADD TO GRAPH FUNCTION #
#########################

function add_to_graphs_blue(graph_list)
    plot!(graph_list[1], (sol.t.-19.0), sol[circ_sys.RV.p], linewidth=2, linecolor=:blue)
    plot!(graph_list[2], (sol.t.-19.0), sol[circ_sys.RA.p], linewidth=2, linecolor=:blue)
    plot!(graph_list[3], (sol.t.-19.0), sol[circ_sys.LV.p], linewidth=2, linecolor=:blue)
    plot!(graph_list[4], (sol.t.-19.0), sol[circ_sys.LA.p], linewidth=2, linecolor=:blue)
    plot!(graph_list[5], sol[circ_sys.RV.V], sol[circ_sys.RV.p], linewidth=2, linecolor=:blue)
    plot!(graph_list[6], sol[circ_sys.RA.V], sol[circ_sys.RA.p], linewidth=2, linecolor=:blue)
    plot!(graph_list[7], sol[circ_sys.LV.V], sol[circ_sys.LV.p], linewidth=2, linecolor=:blue)
    plot!(graph_list[8], sol[circ_sys.LA.V], sol[circ_sys.LA.p], linewidth=2, linecolor=:blue)
    plot!(graph_list[9], (sol.t.-19.0), sol[circ_sys.SAS.C.p], linewidth=2, linecolor=:blue)
    plot!(graph_list[10], (sol.t.-19.0), sol[circ_sys.PAS.C.p], linewidth=2, linecolor=:blue)
    plot!(graph_list[11], (sol.t.-19.0), sol[circ_sys.RV.V], linewidth=2, linecolor=:blue)
    plot!(graph_list[12], (sol.t.-19.0), sol[circ_sys.RA.V], linewidth=2, linecolor=:blue)
    plot!(graph_list[13], (sol.t.-19.0), sol[circ_sys.LV.V], linewidth=2, linecolor=:blue)
    plot!(graph_list[14], (sol.t.-19.0), sol[circ_sys.LA.V], linewidth=2, linecolor=:blue)
    plot!(graph_list[15], (sol.t.-19.0), sol[circ_sys.SAS.C.V], linewidth=2, linecolor=:blue)
    plot!(graph_list[16], (sol.t.-19.0), sol[circ_sys.PAS.C.V], linewidth=2, linecolor=:blue)
end


###################################
# ADD REFERENCE CASE TO ALL PLOTS #
###################################

add_to_graphs_blue([p_rv_pt_ao_r, p_ra_pt_ao_r, p_lv_pt_ao_r, p_la_pt_ao_r, p_rv_pv_ao_r, p_ra_pv_ao_r, p_lv_pv_ao_r, p_la_pv_ao_r,
                p_ao_pt_ao_r, p_pa_pt_ao_r, p_rv_vt_ao_r, p_ra_vt_ao_r, p_lv_vt_ao_r, p_la_vt_ao_r, p_ao_vt_ao_r, p_pa_vt_ao_r]);

add_to_graphs_blue([p_rv_pt_ao_ym, p_ra_pt_ao_ym, p_lv_pt_ao_ym, p_la_pt_ao_ym, p_rv_pv_ao_ym, p_ra_pv_ao_ym, p_lv_pv_ao_ym, p_la_pv_ao_ym,
                p_ao_pt_ao_ym, p_pa_pt_ao_ym, p_rv_vt_ao_ym, p_ra_vt_ao_ym, p_lv_vt_ao_ym, p_la_vt_ao_ym, p_ao_vt_ao_ym, p_pa_vt_ao_ym]);

add_to_graphs_blue([p_rv_pt_pa_r, p_ra_pt_pa_r, p_lv_pt_pa_r, p_la_pt_pa_r, p_rv_pv_pa_r, p_ra_pv_pa_r, p_lv_pv_pa_r, p_la_pv_pa_r,
                p_ao_pt_pa_r, p_pa_pt_pa_r, p_rv_vt_pa_r, p_ra_vt_pa_r, p_lv_vt_pa_r, p_la_vt_pa_r, p_ao_vt_pa_r, p_pa_vt_pa_r]);

add_to_graphs_blue([p_rv_pt_pa_ym, p_ra_pt_pa_ym, p_lv_pt_pa_ym, p_la_pt_pa_ym, p_rv_pv_pa_ym, p_ra_pv_pa_ym, p_lv_pv_pa_ym, p_la_pv_pa_ym,
                p_ao_pt_pa_ym, p_pa_pt_pa_ym, p_rv_vt_pa_ym, p_ra_vt_pa_ym, p_lv_vt_pa_ym, p_la_vt_pa_ym, p_ao_vt_pa_ym, p_pa_vt_pa_ym]);

add_to_graphs_blue([p_rv_pt_la_v0, p_ra_pt_la_v0, p_lv_pt_la_v0, p_la_pt_la_v0, p_rv_pv_la_v0, p_ra_pv_la_v0, p_lv_pv_la_v0, p_la_pv_la_v0,
                p_ao_pt_la_v0, p_pa_pt_la_v0, p_rv_vt_la_v0, p_ra_vt_la_v0, p_lv_vt_la_v0, p_la_vt_la_v0, p_ao_vt_la_v0, p_pa_vt_la_v0]);

add_to_graphs_blue([p_rv_pt_ra_v0_ba, p_ra_pt_ra_v0_ba, p_lv_pt_ra_v0_ba, p_la_pt_ra_v0_ba, p_rv_pv_ra_v0_ba, p_ra_pv_ra_v0_ba, p_lv_pv_ra_v0_ba, p_la_pv_ra_v0_ba,
                p_ao_pt_ra_v0_ba, p_pa_pt_ra_v0_ba, p_rv_vt_ra_v0_ba, p_ra_vt_ra_v0_ba, p_lv_vt_ra_v0_ba, p_la_vt_ra_v0_ba, p_ao_vt_ra_v0_ba, p_pa_vt_ra_v0_ba]);

add_to_graphs_blue([p_rv_pt_lv_emax, p_ra_pt_lv_emax, p_lv_pt_lv_emax, p_la_pt_lv_emax, p_rv_pv_lv_emax, p_ra_pv_lv_emax, p_lv_pv_lv_emax, p_la_pv_lv_emax,
                p_ao_pt_lv_emax, p_pa_pt_lv_emax, p_rv_vt_lv_emax, p_ra_vt_lv_emax, p_lv_vt_lv_emax, p_la_vt_lv_emax, p_ao_vt_lv_emax, p_pa_vt_lv_emax]);

add_to_graphs_blue([p_rv_pt_rv_emax, p_ra_pt_rv_emax, p_lv_pt_rv_emax, p_la_pt_rv_emax, p_rv_pv_rv_emax, p_ra_pv_rv_emax, p_lv_pv_rv_emax, p_la_pv_rv_emax,
                p_ao_pt_rv_emax, p_pa_pt_rv_emax, p_rv_vt_rv_emax, p_ra_vt_rv_emax, p_lv_vt_rv_emax, p_la_vt_rv_emax, p_ao_vt_rv_emax, p_pa_vt_rv_emax]);

add_to_graphs_blue([p_rv_pt_lv_v0, p_ra_pt_lv_v0, p_lv_pt_lv_v0, p_la_pt_lv_v0, p_rv_pv_lv_v0, p_ra_pv_lv_v0, p_lv_pv_lv_v0, p_la_pv_lv_v0,
                p_ao_pt_lv_v0, p_pa_pt_lv_v0, p_rv_vt_lv_v0, p_ra_vt_lv_v0, p_lv_vt_lv_v0, p_la_vt_lv_v0, p_ao_vt_lv_v0, p_pa_vt_lv_v0]);

add_to_graphs_blue([p_rv_pt_ra_v0_bc, p_ra_pt_ra_v0_bc, p_lv_pt_ra_v0_bc, p_la_pt_ra_v0_bc, p_rv_pv_ra_v0_bc, p_ra_pv_ra_v0_bc, p_lv_pv_ra_v0_bc, p_la_pv_ra_v0_bc,
                p_ao_pt_ra_v0_bc, p_pa_pt_ra_v0_bc, p_rv_vt_ra_v0_bc, p_ra_vt_ra_v0_bc, p_lv_vt_ra_v0_bc, p_la_vt_ra_v0_bc, p_ao_vt_ra_v0_bc, p_pa_vt_ra_v0_bc]);

add_to_graphs_blue([p_rv_pt_vc_r, p_ra_pt_vc_r, p_lv_pt_vc_r, p_la_pt_vc_r, p_rv_pv_vc_r, p_ra_pv_vc_r, p_lv_pv_vc_r, p_la_pv_vc_r,
                p_ao_pt_vc_r, p_pa_pt_vc_r, p_rv_vt_vc_r, p_ra_vt_vc_r, p_lv_vt_vc_r, p_la_vt_vc_r, p_ao_vt_vc_r, p_pa_vt_vc_r]);
                
add_to_graphs_blue([p_rv_pt_vc_ym, p_ra_pt_vc_ym, p_lv_pt_vc_ym, p_la_pt_vc_ym, p_rv_pv_vc_ym, p_ra_pv_vc_ym, p_lv_pv_vc_ym, p_la_pv_vc_ym,
                p_ao_pt_vc_ym, p_pa_pt_vc_ym, p_rv_vt_vc_ym, p_ra_vt_vc_ym, p_lv_vt_vc_ym, p_la_vt_vc_ym, p_ao_vt_vc_ym, p_pa_vt_vc_ym]);


#########################
# ADD TO GRAPH FUNCTION #
#########################

function add_to_graphs_red(graph_list)
    plot!(graph_list[1], (sol.t.-19.0), sol[circ_sys.RV.p], linewidth=2, linestyle=:dash, linecolor=:red)
    plot!(graph_list[2], (sol.t.-19.0), sol[circ_sys.RA.p], linewidth=2, linestyle=:dash, linecolor=:red)
    plot!(graph_list[3], (sol.t.-19.0), sol[circ_sys.LV.p], linewidth=2, linestyle=:dash, linecolor=:red)
    plot!(graph_list[4], (sol.t.-19.0), sol[circ_sys.LA.p], linewidth=2, linestyle=:dash, linecolor=:red)
    plot!(graph_list[5], sol[circ_sys.RV.V], sol[circ_sys.RV.p], linewidth=2, linestyle=:dash, linecolor=:red)
    plot!(graph_list[6], sol[circ_sys.RA.V], sol[circ_sys.RA.p], linewidth=2, linestyle=:dash, linecolor=:red)
    plot!(graph_list[7], sol[circ_sys.LV.V], sol[circ_sys.LV.p], linewidth=2, linestyle=:dash, linecolor=:red)
    plot!(graph_list[8], sol[circ_sys.LA.V], sol[circ_sys.LA.p], linewidth=2, linestyle=:dash, linecolor=:red)
    plot!(graph_list[9], (sol.t.-19.0), sol[circ_sys.SAS.C.p], linewidth=2, linestyle=:dash, linecolor=:red)
    plot!(graph_list[10], (sol.t.-19.0), sol[circ_sys.PAS.C.p], linewidth=2, linestyle=:dash, linecolor=:red)
    plot!(graph_list[11], (sol.t.-19.0), sol[circ_sys.RV.V], linewidth=2, linestyle=:dash, linecolor=:red)
    plot!(graph_list[12], (sol.t.-19.0), sol[circ_sys.RA.V], linewidth=2, linestyle=:dash, linecolor=:red)
    plot!(graph_list[13], (sol.t.-19.0), sol[circ_sys.LV.V], linewidth=2, linestyle=:dash, linecolor=:red)
    plot!(graph_list[14], (sol.t.-19.0), sol[circ_sys.LA.V], linewidth=2, linestyle=:dash, linecolor=:red)
    plot!(graph_list[15], (sol.t.-19.0), sol[circ_sys.SAS.C.V], linewidth=2, linestyle=:dash, linecolor=:red)
    plot!(graph_list[16], (sol.t.-19.0), sol[circ_sys.PAS.C.V], linewidth=2, linestyle=:dash, linecolor=:red)
end


###################################
# DECREASE IN AORTIC SINUS RADIUS #
###################################

# Change necessary parameters
@named SAS = CRL(C=Csas*0.125, R=Rsas*16, L=Lsas*4);

# Create ODE system
@named _circ_model = ODESystem(circ_eqs, t)
@named circ_model = compose(_circ_model,
[LV, RV, LA, RA, AV, MV, PV, TV, SAS, SAT, SAR, SCP, SVN, PAS, PAT, PAR, PCP, PVN])
circ_sys = structural_simplify(circ_model)

# Solve system of ODEs
prob = ODAEProblem(circ_sys, u0, (0.0, 20.0))
sol = solve(prob, RK4(), reltol=1e-9, abstol=1e-12, saveat=19:0.01:20)

# Add to graphs
add_to_graphs_red([p_rv_pt_ao_r, p_ra_pt_ao_r, p_lv_pt_ao_r, p_la_pt_ao_r, p_rv_pv_ao_r, p_ra_pv_ao_r, p_lv_pv_ao_r, p_la_pv_ao_r,
                p_ao_pt_ao_r, p_pa_pt_ao_r, p_rv_vt_ao_r, p_ra_vt_ao_r, p_lv_vt_ao_r, p_la_vt_ao_r, p_ao_vt_ao_r, p_pa_vt_ao_r]);

# Revert to default case
@named SAS = CRL(C=Csas, R=Rsas, L=Lsas);


############################################
# DECREASE IN AORTIC SINUS YOUNG'S MODULUS #
############################################

# Change necessary parameters
@named SAS = CRL(C=Csas*3, R=Rsas, L=Lsas);

# Create ODE system
@named _circ_model = ODESystem(circ_eqs, t)
@named circ_model = compose(_circ_model,
[LV, RV, LA, RA, AV, MV, PV, TV, SAS, SAT, SAR, SCP, SVN, PAS, PAT, PAR, PCP, PVN])
circ_sys = structural_simplify(circ_model)

# Solve system of ODEs
prob = ODAEProblem(circ_sys, u0, (0.0, 20.0))
sol = solve(prob, RK4(), reltol=1e-9, abstol=1e-12, saveat=19:0.01:20)

# Add to graphs
add_to_graphs_red([p_rv_pt_ao_ym, p_ra_pt_ao_ym, p_lv_pt_ao_ym, p_la_pt_ao_ym, p_rv_pv_ao_ym, p_ra_pv_ao_ym, p_lv_pv_ao_ym, p_la_pv_ao_ym,
                p_ao_pt_ao_ym, p_pa_pt_ao_ym, p_rv_vt_ao_ym, p_ra_vt_ao_ym, p_lv_vt_ao_ym, p_la_vt_ao_ym, p_ao_vt_ao_ym, p_pa_vt_ao_ym]);

# Revert to default case
@named SAS = CRL(C=Csas, R=Rsas, L=Lsas);


######################################
# DECREASE IN PULMONARY SINUS RADIUS #
######################################

# Change necessary parameters
@named PAS = CRL(C=Cpas*0.204, R=Rpas*8.352, L=Lpas*2.89);

# Create ODE system
@named _circ_model = ODESystem(circ_eqs, t)
@named circ_model = compose(_circ_model,
[LV, RV, LA, RA, AV, MV, PV, TV, SAS, SAT, SAR, SCP, SVN, PAS, PAT, PAR, PCP, PVN])
circ_sys = structural_simplify(circ_model)

# Solve system of ODEs
prob = ODAEProblem(circ_sys, u0, (0.0, 20.0))
sol = solve(prob, RK4(), reltol=1e-9, abstol=1e-12, saveat=19:0.01:20)

# Add to graphs
add_to_graphs_red([p_rv_pt_pa_r, p_ra_pt_pa_r, p_lv_pt_pa_r, p_la_pt_pa_r, p_rv_pv_pa_r, p_ra_pv_pa_r, p_lv_pv_pa_r, p_la_pv_pa_r,
                p_ao_pt_pa_r, p_pa_pt_pa_r, p_rv_vt_pa_r, p_ra_vt_pa_r, p_lv_vt_pa_r, p_la_vt_pa_r, p_ao_vt_pa_r, p_pa_vt_pa_r]);

# Revert to default
@named PAS = CRL(C=Cpas, R=Rpas, L=Lpas);


#############################################
# DECREASED PULMONARY SINUS YOUNG'S MODULUS #
#############################################

# Change necessary parameters
@named PAS = CRL(C=Cpas*3, R=Rpas, L=Lpas);

# Create ODE system
@named _circ_model = ODESystem(circ_eqs, t)
@named circ_model = compose(_circ_model,
[LV, RV, LA, RA, AV, MV, PV, TV, SAS, SAT, SAR, SCP, SVN, PAS, PAT, PAR, PCP, PVN])
circ_sys = structural_simplify(circ_model)

# Solve system of ODEs
prob = ODAEProblem(circ_sys, u0, (0.0, 20.0))
sol = solve(prob, RK4(), reltol=1e-9, abstol=1e-12, saveat=19:0.01:20)

# Add to graphs
add_to_graphs_red([p_rv_pt_pa_ym, p_ra_pt_pa_ym, p_lv_pt_pa_ym, p_la_pt_pa_ym, p_rv_pv_pa_ym, p_ra_pv_pa_ym, p_lv_pv_pa_ym, p_la_pv_pa_ym,
                p_ao_pt_pa_ym, p_pa_pt_pa_ym, p_rv_vt_pa_ym, p_ra_vt_pa_ym, p_lv_vt_pa_ym, p_la_vt_pa_ym, p_ao_vt_pa_ym, p_pa_vt_pa_ym]);

# Revert to default
@named PAS = CRL(C=Cpas, R=Rpas, L=Lpas);


################################
# DECREASED LEFT ATRIUM VOLUME #
################################

# Change necessary parameters
@named LA = HeartChamber(V_0=v0_la*0.75, p_0=p0_la, E_min=Emin_la, E_max=Emax_la, T=T, T_es=Tpww_la / 2, T_ep=Tpww_la, Eshift=Tpwb_la);

# Create ODE system
@named _circ_model = ODESystem(circ_eqs, t)
@named circ_model = compose(_circ_model,
[LV, RV, LA, RA, AV, MV, PV, TV, SAS, SAT, SAR, SCP, SVN, PAS, PAT, PAR, PCP, PVN])
circ_sys = structural_simplify(circ_model)

# Solve system of ODEs
prob = ODAEProblem(circ_sys, u0, (0.0, 20.0))
sol = solve(prob, RK4(), reltol=1e-9, abstol=1e-12, saveat=19:0.01:20)

# Add to graphs
add_to_graphs_red([p_rv_pt_la_v0, p_ra_pt_la_v0, p_lv_pt_la_v0, p_la_pt_la_v0, p_rv_pv_la_v0, p_ra_pv_la_v0, p_lv_pv_la_v0, p_la_pv_la_v0,
                p_ao_pt_la_v0, p_pa_pt_la_v0, p_rv_vt_la_v0, p_ra_vt_la_v0, p_lv_vt_la_v0, p_la_vt_la_v0, p_ao_vt_la_v0, p_pa_vt_la_v0]);

# Revert to default
@named LA = HeartChamber(V_0=v0_la, p_0=p0_la, E_min=Emin_la, E_max=Emax_la, T=T, T_es=Tpww_la / 2, T_ep=Tpww_la, Eshift=Tpwb_la);


###########################################
# DECREASED LEFT ATRIUM VOLUME (BIATRIAL) #
###########################################

# Change necessary parameters
@named RA = HeartChamber(V_0=v0_ra*0.75, p_0=p0_ra, E_min=Emin_ra, E_max=Emax_ra, T=T, T_es=Tpww_ra / 2, T_ep=Tpww_ra, Eshift=Tpwb_ra);

# Create ODE system
@named _circ_model = ODESystem(circ_eqs, t)
@named circ_model = compose(_circ_model,
[LV, RV, LA, RA, AV, MV, PV, TV, SAS, SAT, SAR, SCP, SVN, PAS, PAT, PAR, PCP, PVN])
circ_sys = structural_simplify(circ_model)

# Solve system of ODEs
prob = ODAEProblem(circ_sys, u0, (0.0, 20.0))
sol = solve(prob, RK4(), reltol=1e-9, abstol=1e-12, saveat=19:0.01:20)

# Add to graphs
add_to_graphs_red([p_rv_pt_ra_v0_ba, p_ra_pt_ra_v0_ba, p_lv_pt_ra_v0_ba, p_la_pt_ra_v0_ba, p_rv_pv_ra_v0_ba, p_ra_pv_ra_v0_ba, p_lv_pv_ra_v0_ba, p_la_pv_ra_v0_ba,
                p_ao_pt_ra_v0_ba, p_pa_pt_ra_v0_ba, p_rv_vt_ra_v0_ba, p_ra_vt_ra_v0_ba, p_lv_vt_ra_v0_ba, p_la_vt_ra_v0_ba, p_ao_vt_ra_v0_ba, p_pa_vt_ra_v0_ba]);

# Revert to default
@named RA = HeartChamber(V_0=v0_ra, p_0=p0_ra, E_min=Emin_ra, E_max=Emax_ra, T=T, T_es=Tpww_ra / 2, T_ep=Tpww_ra, Eshift=Tpwb_ra);


######################################
# INCREASED LEFT VENTRICLE ELASTANCE #
######################################

# cHANGE NECESSARY PARAMETERS
@named LV = HeartChamber(V_0=v0_lv, p_0=p0_lv, E_min=Emin_lv, E_max=Emax_lv*1.75, T=T, T_es=Tes_lv, T_ep=Ted_lv, Eshift=0.0);

# Create ODE system
@named _circ_model = ODESystem(circ_eqs, t)
@named circ_model = compose(_circ_model,
[LV, RV, LA, RA, AV, MV, PV, TV, SAS, SAT, SAR, SCP, SVN, PAS, PAT, PAR, PCP, PVN])
circ_sys = structural_simplify(circ_model)

# Solve system of ODEs
prob = ODAEProblem(circ_sys, u0, (0.0, 20.0))
sol = solve(prob, RK4(), reltol=1e-9, abstol=1e-12, saveat=19:0.01:20)

# Add to plots
add_to_graphs_red([p_rv_pt_lv_emax, p_ra_pt_lv_emax, p_lv_pt_lv_emax, p_la_pt_lv_emax, p_rv_pv_lv_emax, p_ra_pv_lv_emax, p_lv_pv_lv_emax, p_la_pv_lv_emax,
                p_ao_pt_lv_emax, p_pa_pt_lv_emax, p_rv_vt_lv_emax, p_ra_vt_lv_emax, p_lv_vt_lv_emax, p_la_vt_lv_emax, p_ao_vt_lv_emax, p_pa_vt_lv_emax]);

#Revert to default
@named LV = HeartChamber(V_0=v0_lv, p_0=p0_lv, E_min=Emin_lv, E_max=Emax_lv, T=T, T_es=Tes_lv, T_ep=Ted_lv, Eshift=0.0);


#######################################
# INCREASED RIGHT VENTRICLE ELASTANCE #
#######################################

# Change necessary parameters
@named RV = HeartChamber(V_0=v0_rv, p_0=p0_rv, E_min=Emin_rv, E_max=Emax_rv*1.1, T=T, T_es=Tes_rv, T_ep=Ted_rv, Eshift=0.0);

# Create ODE system
@named _circ_model = ODESystem(circ_eqs, t)
@named circ_model = compose(_circ_model,
[LV, RV, LA, RA, AV, MV, PV, TV, SAS, SAT, SAR, SCP, SVN, PAS, PAT, PAR, PCP, PVN])
circ_sys = structural_simplify(circ_model)

# Solve system of ODEs
prob = ODAEProblem(circ_sys, u0, (0.0, 20.0))
sol = solve(prob, RK4(), reltol=1e-9, abstol=1e-12, saveat=19:0.01:20)

# Add to graphs 
add_to_graphs_red([p_rv_pt_rv_emax, p_ra_pt_rv_emax, p_lv_pt_rv_emax, p_la_pt_rv_emax, p_rv_pv_rv_emax, p_ra_pv_rv_emax, p_lv_pv_rv_emax, p_la_pv_rv_emax,
                p_ao_pt_rv_emax, p_pa_pt_rv_emax, p_rv_vt_rv_emax, p_ra_vt_rv_emax, p_lv_vt_rv_emax, p_la_vt_rv_emax, p_ao_vt_rv_emax, p_pa_vt_rv_emax]);

# Revert to default
@named RV = HeartChamber(V_0=v0_rv, p_0=p0_rv, E_min=Emin_rv, E_max=Emax_rv, T=T, T_es=Tes_rv, T_ep=Ted_rv, Eshift=0.0);


###################################
# INCREASED LEFT VENTRICLE VOLUME #
###################################

# Change necessary parameters
@named LV = HeartChamber(V_0=v0_lv*1.5, p_0=p0_lv, E_min=Emin_lv, E_max=Emax_lv, T=T, T_es=Tes_lv, T_ep=Ted_lv, Eshift=0.0);

# Create ODE system
@named _circ_model = ODESystem(circ_eqs, t)
@named circ_model = compose(_circ_model,
[LV, RV, LA, RA, AV, MV, PV, TV, SAS, SAT, SAR, SCP, SVN, PAS, PAT, PAR, PCP, PVN])
circ_sys = structural_simplify(circ_model)

# Solve system of ODEs
prob = ODAEProblem(circ_sys, u0, (0.0, 20.0))
sol = solve(prob, RK4(), reltol=1e-9, abstol=1e-12, saveat=19:0.01:20)

# Add to graphs 
add_to_graphs_red([p_rv_pt_lv_v0, p_ra_pt_lv_v0, p_lv_pt_lv_v0, p_la_pt_lv_v0, p_rv_pv_lv_v0, p_ra_pv_lv_v0, p_lv_pv_lv_v0, p_la_pv_lv_v0,
                p_ao_pt_lv_v0, p_pa_pt_lv_v0, p_rv_vt_lv_v0, p_ra_vt_lv_v0, p_lv_vt_lv_v0, p_la_vt_lv_v0, p_ao_vt_lv_v0, p_pa_vt_lv_v0]);

# Revert to default
@named LV = HeartChamber(V_0=v0_lv, p_0=p0_lv, E_min=Emin_lv, E_max=Emax_lv, T=T, T_es=Tes_lv, T_ep=Ted_lv, Eshift=0.0);


###########################################
# DECREASED RIGHT ATRIUM VOLUME (BICAVAL) #
###########################################

# Change necessary parameters
@named RA = HeartChamber(V_0=v0_ra*0.5, p_0=p0_ra, E_min=Emin_ra, E_max=Emax_ra, T=T, T_es=Tpww_ra / 2, T_ep=Tpww_ra, Eshift=Tpwb_ra);

# Create ODE system
@named _circ_model = ODESystem(circ_eqs, t)
@named circ_model = compose(_circ_model,
[LV, RV, LA, RA, AV, MV, PV, TV, SAS, SAT, SAR, SCP, SVN, PAS, PAT, PAR, PCP, PVN])
circ_sys = structural_simplify(circ_model)

# Solve system of ODEs
prob = ODAEProblem(circ_sys, u0, (0.0, 20.0))
sol = solve(prob, RK4(), reltol=1e-9, abstol=1e-12, saveat=19:0.01:20)

# Add to graphs 
add_to_graphs_red([p_rv_pt_ra_v0_bc, p_ra_pt_ra_v0_bc, p_lv_pt_ra_v0_bc, p_la_pt_ra_v0_bc, p_rv_pv_ra_v0_bc, p_ra_pv_ra_v0_bc, p_lv_pv_ra_v0_bc, p_la_pv_ra_v0_bc,
                p_ao_pt_ra_v0_bc, p_pa_pt_ra_v0_bc, p_rv_vt_ra_v0_bc, p_ra_vt_ra_v0_bc, p_lv_vt_ra_v0_bc, p_la_vt_ra_v0_bc, p_ao_vt_ra_v0_bc, p_pa_vt_ra_v0_bc]);

# Revert to default
@named RA = HeartChamber(V_0=v0_ra, p_0=p0_ra, E_min=Emin_ra, E_max=Emax_ra, T=T, T_es=Tpww_ra / 2, T_ep=Tpww_ra, Eshift=Tpwb_ra);


##############################
# DECREASED VENA CAVA RADIUS #
##############################

# Change necessary parameters
@named SVN = CR(R=(0.075*0.995)+(0.075*0.005*1.916), C=(20.5*0.95)+(20.5*0.05*0.614));

# Create ODE system
@named _circ_model = ODESystem(circ_eqs, t)
@named circ_model = compose(_circ_model,
[LV, RV, LA, RA, AV, MV, PV, TV, SAS, SAT, SAR, SCP, SVN, PAS, PAT, PAR, PCP, PVN])
circ_sys = structural_simplify(circ_model)

# Solve system of ODEs
prob = ODAEProblem(circ_sys, u0, (0.0, 20.0))
sol = solve(prob, RK4(), reltol=1e-9, abstol=1e-12, saveat=19:0.01:20)

# Add to graphs
add_to_graphs_red([p_rv_pt_vc_r, p_ra_pt_vc_r, p_lv_pt_vc_r, p_la_pt_vc_r, p_rv_pv_vc_r, p_ra_pv_vc_r, p_lv_pv_vc_r, p_la_pv_vc_r,
                p_ao_pt_vc_r, p_pa_pt_vc_r, p_rv_vt_vc_r, p_ra_vt_vc_r, p_lv_vt_vc_r, p_la_vt_vc_r, p_ao_vt_vc_r, p_pa_vt_vc_r]);

# Revert to default
@named SVN = CR(R=Rsvn, C=Csvn);


#######################################
# DECREASED VENA CAVA YOUNG'S MODULUS #
#######################################

# Change necessary parameters
@named SVN = CR(R=Rsvn, C=(20.5*0.95)+(20.5*0.05*3));

# Create ODE system
@named _circ_model = ODESystem(circ_eqs, t)
@named circ_model = compose(_circ_model,
[LV, RV, LA, RA, AV, MV, PV, TV, SAS, SAT, SAR, SCP, SVN, PAS, PAT, PAR, PCP, PVN])
circ_sys = structural_simplify(circ_model)

# Solve system of ODEs
prob = ODAEProblem(circ_sys, u0, (0.0, 20.0))
sol = solve(prob, RK4(), reltol=1e-9, abstol=1e-12, saveat=19:0.01:20)

# Add to graphs 
add_to_graphs_red([p_rv_pt_vc_ym, p_ra_pt_vc_ym, p_lv_pt_vc_ym, p_la_pt_vc_ym, p_rv_pv_vc_ym, p_ra_pv_vc_ym, p_lv_pv_vc_ym, p_la_pv_vc_ym,
                p_ao_pt_vc_ym, p_pa_pt_vc_ym, p_rv_vt_vc_ym, p_ra_vt_vc_ym, p_lv_vt_vc_ym, p_la_vt_vc_ym, p_ao_vt_vc_ym, p_pa_vt_vc_ym]);

# Revert to default
@named SVN = CR(R=Rsvn, C=Csvn);


#############################
# MAKING AND SAVING FIGURES #
#############################

# Pressure in Left Ventricle
fig = plot(p_lv_pt_ao_r, p_lv_pt_ao_ym, p_lv_pt_pa_r, 
        p_lv_pt_pa_ym, p_lv_pt_la_v0, p_lv_pt_ra_v0_ba, 
        p_lv_pt_lv_emax, p_lv_pt_rv_emax, p_lv_pt_lv_v0,
        p_lv_pt_ra_v0_bc, p_lv_pt_vc_r, p_lv_pt_vc_ym, 
        layout=(4,3), margin=10Plots.mm)
plot!(size=(1200,1200))
savefig(fig, "plots/lv_pt.pdf")

# Pressure in Right Ventricle
fig = plot(p_rv_pt_ao_r, p_rv_pt_ao_ym, p_rv_pt_pa_r, 
        p_rv_pt_pa_ym, p_rv_pt_la_v0, p_rv_pt_ra_v0_ba, 
        p_rv_pt_lv_emax, p_rv_pt_rv_emax, p_rv_pt_lv_v0,
        p_rv_pt_ra_v0_bc, p_rv_pt_vc_r, p_rv_pt_vc_ym, 
        layout=(4,3), margin=10Plots.mm)
plot!(size=(1200,1200))
savefig(fig, "plots/rv_pt.pdf")

# Pressure in Left Atrium
fig = plot(p_la_pt_ao_r, p_la_pt_ao_ym, p_la_pt_pa_r, 
        p_la_pt_pa_ym, p_la_pt_la_v0, p_la_pt_ra_v0_ba, 
        p_la_pt_lv_emax, p_la_pt_rv_emax, p_la_pt_lv_v0,
        p_la_pt_ra_v0_bc, p_la_pt_vc_r, p_la_pt_vc_ym, 
        layout=(4,3), margin=10Plots.mm)
plot!(size=(1200,1200))
savefig(fig, "plots/la_pt.pdf")

# Pressure in right Atrium
fig = plot(p_ra_pt_ao_r, p_ra_pt_ao_ym, p_ra_pt_pa_r, 
        p_ra_pt_pa_ym, p_ra_pt_la_v0, p_ra_pt_ra_v0_ba, 
        p_ra_pt_lv_emax, p_ra_pt_rv_emax, p_ra_pt_lv_v0,
        p_ra_pt_ra_v0_bc, p_ra_pt_vc_r, p_ra_pt_vc_ym, 
        layout=(4,3), margin=10Plots.mm)
plot!(size=(1200,1200))
savefig(fig, "plots/ra_pt.pdf")

# Left Ventricle PV loop
fig = plot(p_lv_pv_ao_r, p_lv_pv_ao_ym, p_lv_pv_pa_r, 
        p_lv_pv_pa_ym, p_lv_pv_la_v0, p_lv_pv_ra_v0_ba, 
        p_lv_pv_lv_emax, p_lv_pv_rv_emax, p_lv_pv_lv_v0,
        p_lv_pv_ra_v0_bc, p_lv_pv_vc_r, p_lv_pv_vc_ym, 
        layout=(4,3), margin=10Plots.mm)
plot!(size=(1200,1200))
savefig(fig, "plots/lv_pv.pdf")

# Right Ventricle PV loop
fig = plot(p_rv_pv_ao_r, p_rv_pv_ao_ym, p_rv_pv_pa_r, 
        p_rv_pv_pa_ym, p_rv_pv_la_v0, p_rv_pv_ra_v0_ba, 
        p_rv_pv_lv_emax, p_rv_pv_rv_emax, p_rv_pv_lv_v0,
        p_rv_pv_ra_v0_bc, p_rv_pv_vc_r, p_rv_pv_vc_ym, 
        layout=(4,3), margin=10Plots.mm)
plot!(size=(1200,1200))
savefig(fig, "plots/rv_pv.pdf")

# Left Atrium PV loop
fig = plot(p_la_pv_ao_r, p_la_pv_ao_ym, p_la_pv_pa_r, 
        p_la_pv_pa_ym, p_la_pv_la_v0, p_la_pv_ra_v0_ba, 
        p_la_pv_lv_emax, p_la_pv_rv_emax, p_la_pv_lv_v0,
        p_la_pv_ra_v0_ba, p_la_pv_vc_r, p_la_pv_vc_ym, 
        layout=(4,3), margin=10Plots.mm)
plot!(size=(1200,1200))
savefig(fig, "plots/la_pv.pdf")

# Right Atrium PV loop
fig = plot(p_ra_pv_ao_r, p_ra_pv_ao_ym, p_ra_pv_pa_r, 
        p_ra_pv_pa_ym, p_ra_pv_la_v0, p_ra_pv_ra_v0_ba, 
        p_ra_pv_lv_emax, p_ra_pv_rv_emax, p_ra_pv_lv_v0,
        p_ra_pv_ra_v0_bc, p_ra_pv_vc_r, p_ra_pv_vc_ym, 
        layout=(4,3), margin=10Plots.mm)
plot!(size=(1200,1200))
savefig(fig, "plots/ra_pv.pdf")

# Pressure in Aortic sinus
fig = plot(p_ao_pt_ao_r, p_ao_pt_ao_ym, p_ao_pt_pa_r, 
        p_ao_pt_pa_ym, p_ao_pt_la_v0, p_ao_pt_ra_v0_ba, 
        p_ao_pt_lv_emax, p_ao_pt_rv_emax, p_ao_pt_lv_v0,
        p_ao_pt_ra_v0_bc, p_ao_pt_vc_r, p_ao_pt_vc_ym, 
        layout=(4,3), margin=10Plots.mm)
plot!(size=(1200,1200))
savefig(fig, "plots/ao_pt.pdf")

# Pressure in Pulmonary Sinus
fig = plot(p_pa_pt_ao_r, p_ao_pt_ao_ym, p_pa_pt_pa_r, 
        p_pa_pt_pa_ym, p_pa_pt_la_v0, p_pa_pt_ra_v0_ba, 
        p_pa_pt_lv_emax, p_pa_pt_rv_emax, p_pa_pt_lv_v0,
        p_pa_pt_ra_v0_bc, p_pa_pt_vc_r, p_pa_pt_vc_ym, 
        layout=(4,3), margin=10Plots.mm)
plot!(size=(1200,1200))
savefig(fig, "plots/pa_pt.pdf")

# Volume in Left Ventricle
fig = plot(p_lv_vt_ao_r, p_lv_vt_ao_ym, p_lv_vt_pa_r, 
        p_lv_vt_pa_ym, p_lv_vt_la_v0, p_lv_vt_ra_v0_ba, 
        p_lv_vt_lv_emax, p_lv_vt_rv_emax, p_lv_vt_lv_v0,
        p_lv_vt_ra_v0_bc, p_lv_vt_vc_r, p_lv_vt_vc_ym, 
        layout=(4,3), margin=10Plots.mm)
plot!(size=(1200,1200))
savefig(fig, "plots/lv_vt.pdf")

# Volume in Right Ventricle
fig = plot(p_rv_vt_ao_r, p_rv_vt_ao_ym, p_rv_vt_pa_r, 
        p_rv_vt_pa_ym, p_rv_vt_la_v0, p_rv_vt_ra_v0_ba, 
        p_rv_vt_lv_emax, p_rv_vt_rv_emax, p_rv_vt_lv_v0,
        p_rv_vt_ra_v0_bc, p_rv_vt_vc_r, p_rv_vt_vc_ym, 
        layout=(4,3), margin=10Plots.mm)
plot!(size=(1200,1200))
savefig(fig, "plots/rv_vt.pdf")

# Volume in Left Atrium
fig = plot(p_la_vt_ao_r, p_la_vt_ao_ym, p_la_vt_pa_r, 
        p_la_vt_pa_ym, p_la_vt_la_v0, p_la_vt_ra_v0_ba, 
        p_la_vt_lv_emax, p_la_vt_rv_emax, p_la_vt_lv_v0,
        p_la_vt_ra_v0_bc, p_la_vt_vc_r, p_la_vt_vc_ym, 
        layout=(4,3), margin=10Plots.mm)
plot!(size=(1200,1200))
savefig(fig, "plots/la_vt.pdf")

# Volume in Right Atrium
fig = plot(p_ra_vt_ao_r, p_ra_vt_ao_ym, p_ra_vt_pa_r, 
        p_ra_vt_pa_ym, p_ra_vt_la_v0, p_ra_vt_ra_v0_ba, 
        p_ra_vt_lv_emax, p_ra_vt_rv_emax, p_ra_vt_lv_v0,
        p_ra_vt_ra_v0_bc, p_ra_vt_vc_r, p_ra_vt_vc_ym, 
        layout=(4,3), margin=10Plots.mm)
plot!(size=(1200,1200))
savefig(fig, "plots/ra_vt.pdf")

# Volume in Aortic Sinus
fig = plot(p_ao_vt_ao_r, p_ao_vt_ao_ym, p_ao_vt_pa_r, 
        p_ao_vt_pa_ym, p_ao_vt_la_v0, p_ao_vt_ra_v0_ba, 
        p_ao_vt_lv_emax, p_ao_vt_rv_emax, p_ao_vt_lv_v0,
        p_ao_vt_ra_v0_bc, p_ao_vt_vc_r, p_ao_vt_vc_ym, 
        layout=(4,3), margin=10Plots.mm)
plot!(size=(1200,1200))
savefig(fig, "plots/ao_vt.pdf")

# Volume in Pulmonary Sinus
fig = plot(p_pa_vt_ao_r, p_pa_vt_ao_ym, p_pa_vt_pa_r, 
        p_pa_vt_pa_ym, p_pa_vt_la_v0, p_pa_vt_ra_v0_ba, 
        p_pa_vt_lv_emax, p_pa_vt_rv_emax, p_pa_vt_lv_v0,
        p_pa_vt_ra_v0_bc, p_pa_vt_vc_r, p_pa_vt_vc_ym,
        layout=(4,3), margin=10Plots.mm)
plot!(size=(1200,1200))
savefig(fig, "plots/pa_vt.pdf")