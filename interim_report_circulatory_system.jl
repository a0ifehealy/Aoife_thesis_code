using ModelingToolkit
using DifferentialEquations
using Plots

@variables t
D = Differential(t)
include("model_parameters.jl")
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


###########################################
# DEFINE COMPONENTS OF CIRCULATORY SYSTEM #
###########################################

# Heart chambers
@named LV = HeartChamber(V₀=v0_lv, p₀=p0_lv, E_min=Emin_lv, E_max=Emax_lv, T=T, T_es=Tes_lv, T_ep=Ted_lv, Eshift=0.0)
@named LA = HeartChamber(V₀=v0_la, p₀=p0_la, E_min=Emin_la, E_max=Emax_la, T=T, T_es=Tpww_la / 2, T_ep=Tpww_la, Eshift=Tpwb_la)
@named RV = HeartChamber(V₀=v0_rv, p₀=p0_rv, E_min=Emin_rv, E_max=Emax_rv, T=T, T_es=Tes_rv, T_ep=Ted_rv, Eshift=0.0)
@named RA = HeartChamber(V₀=v0_ra, p₀=p0_ra, E_min=Emin_ra, E_max=Emax_ra, T=T, T_es=Tpww_ra / 2, T_ep=Tpww_ra, Eshift=Tpwb_ra)

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
@named PVN = CR(R=Rpvn, C=Cpvn)


################################
# COMPOSE AND SOLVE ODE SYSTEM #
################################

# Connect components to form circuit
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
sol = solve(prob, RK4(), reltol=1e-9, abstol=1e-12, saveat=19:0.01:20)
plot(sol)


####################
# PLOTTING RESULTS #
####################

## ELASTANCE OF CHAMBERS OVER CYCLE ##

t1 = 0:0.01:1
e_lv = []; e_rv = []; e_la = []; e_ra = []
for t in t1
    append!(e_lv, Elastance(t, Emin_lv, Emax_lv, T, Tes_lv, Ted_lv, Eshift_lv))
    append!(e_rv, Elastance(t, Emin_rv, Emax_rv, T, Tes_rv, Ted_rv, Eshift_rv))
    append!(e_la, Elastance(t, Emin_la, Emax_la, T, Tes_la, Ted_la, Eshift_la))
    append!(e_ra, Elastance(t, Emin_ra, Emax_ra, T, Tes_ra, Ted_ra, Eshift_ra))
end

fig = plot(t1, e_lv, labels="Left Ventricle", ylabel="Elastance (mmHg/ml)", yguidefontsize=20, ytickfontsize=16, 
            xlabel="t (s)", xguidefontsize=20, xtickfontsize=16, legendfontsize=16, lw=2, margin=5Plots.mm)
plot!(t1, e_rv, labels="Right Ventricle", ls=:dash, lw=2)
plot!(t1, e_la, labels="Left Atrium", ls=:dot, lw=2)
plot!(t1, e_ra, labels="Right Atrium", ls=:dashdot, lw=2)
plot!(size=(1200,700))
savefig(fig, "plots/elastance")


## RE-INDEXING RESULTS SO THAT THE TIME STARTS AFTER THE SOLUTION HAS STABILISED ##

t = sol.t .- sol.t[1];


## PRESSURE AND ELASTANCE OF HEART CHAMBERS ##

# Left Ventricle
p_lv = plot(t, sol[circ_sys.LV.p], c=:blue, xlabel="t (s)", ylabel="Pressure (mmHg)", yguidefontsize=12, 
            ytickfontsize=12, xguidefontsize=12, xtickfontsize=12, framestyle=:zerolines, title = "Left Ventricle", 
            margin=5Plots.mm, label = "")
plot!(twinx(), t, e_lv, c=:gray30, ls=:dash, grid=true, xlabel="", ylabel="Elastance (mmHg/mL)", label="",
        yguidefontsize=12, ytickfontsize=12);

# Right Ventricle
p_rv = plot(t, sol[circ_sys.RV.p], c=:blue, xlabel="t (s)", ylabel="Pressure (mmHg)", yguidefontsize=12, 
            ytickfontsize=12, xguidefontsize=12, xtickfontsize=12, framestyle=:zerolines, title = "Right Ventricle", 
            margin=5Plots.mm, label = "")
plot!(twinx(), t, e_rv, c=:gray30, ls=:dash, grid=true, xlabel="", ylabel="Elastance (mmHg/mL)", label="",
        yguidefontsize=12, ytickfontsize=12);

# Left Atrium
p_la = plot(t, sol[circ_sys.LA.p], c=:blue, xlabel="t (s)", ylabel="Pressure (mmHg)", yguidefontsize=12, 
            ytickfontsize=12, xguidefontsize=12, xtickfontsize=12, framestyle=:zerolines, title = "Left Atrium", 
            margin=5Plots.mm, label = "")
plot!(twinx(), t, e_la, c=:gray30, ls=:dash, grid=true, xlabel="", ylabel="Elastance (mmHg/mL)", label="",  
        yguidefontsize=12, ytickfontsize=12);
        
# Right Atrium
p_ra = plot(t, sol[circ_sys.RA.p], c=:blue, xlabel="t (s)", ylabel="Pressure (mmHg)", yguidefontsize=12, 
            ytickfontsize=12, xguidefontsize=12, xtickfontsize=12, framestyle=:zerolines, title = "Right Atrium", 
            margin=5Plots.mm, label = "")
plot!(twinx(), t, e_ra, c=:gray30, ls=:dash, grid=true, xlabel="", ylabel="Elastance (mmHg/mL)", label="",  
        yguidefontsize=12, ytickfontsize=12);

# Adding four plots to one figure
fig = plot(p_lv, p_rv, p_la, p_ra, layout=(2,2), margin=5Plots.mm, legend = false)
plot!(size=(1200,700))
savefig(fig, "plots/chambers_pressure_elastance")

# Create legend
legend = plot([0], showaxis = false, grid = false, label = "Pressure in chamber", legend = :topleft, margin=5Plots.mm, 
            legendfontsize=16, c=:blue)
plot!([0], showaxis = false, grid = false, label = "Elastance of chamber walls", c=:gray30, ls=:dash);   
savefig(legend, "plots/chambers_pressure_elastance_legend") 


## VOLUME AND ELASTANCE OF HEART CHAMBERS ##

# Left Ventricle
p_lv = plot(t, sol[circ_sys.LV.V], c=:red, xlabel="t (s)", ylabel="Volume (ml)", yguidefontsize=12, ytickfontsize=12,
            xguidefontsize=12, xtickfontsize=12, framestyle=:zerolines, title = "Left Ventricle", margin=5Plots.mm, 
            label = "")
plot!(twinx(), t, e_lv, c=:gray30, ls=:dash, grid=true, xlabel="", ylabel="Elastance (mmHg/mL)", label="",  
        yguidefontsize=12, ytickfontsize=12);

# Right Ventricle
p_rv = plot(t, sol[circ_sys.RV.V], c=:red, xlabel="t (s)", ylabel="Volume (ml)", yguidefontsize=12, ytickfontsize=12, 
            xguidefontsize=12, xtickfontsize=12, framestyle=:zerolines, title = "Right Ventricle", margin=5Plots.mm, 
            label = "")
plot!(twinx(), t, e_rv, c=:gray30, ls=:dash, grid=true, xlabel="", ylabel="Elastance (mmHg/mL)", label="",  
        yguidefontsize=12, ytickfontsize=12);

# Left Atrium
p_la = plot(t, sol[circ_sys.LA.V], c=:red, xlabel="t (s)", ylabel="Volume (ml)", yguidefontsize=12, ytickfontsize=12, 
            xguidefontsize=12, xtickfontsize=12, framestyle=:zerolines, title = "Left Atrium", margin=5Plots.mm, 
            label = "")
plot!(twinx(), t, e_la, c=:gray30, ls=:dash, grid=true, xlabel="", ylabel="Elastance (mmHg/mL)", label="",  
        yguidefontsize=12, ytickfontsize=12);

# Right Atrium       
p_ra = plot(t, sol[circ_sys.RA.V], c=:red, xlabel="t (s)", ylabel="Volume (ml)", yguidefontsize=12, ytickfontsize=12, 
            xguidefontsize=12, xtickfontsize=12, framestyle=:zerolines, title = "Right Atrium", margin=5Plots.mm, 
            label = "")
plot!(twinx(), t, e_ra, c=:gray30, ls=:dash, grid=true, xlabel="", ylabel="Elastance (mmHg/mL)", label="",  
        yguidefontsize=12, ytickfontsize=12);

# Adding four plots to one figure
fig = plot(p_lv, p_rv, p_la, p_ra, layout=(2,2), margin=5Plots.mm, legend = false)
plot!(size=(1200,700))
savefig(fig, "plots/chambers_volume_elastance")

# Create legend
legend = plot([0], showaxis = false, grid = false, label = "Volume of blood in chamber", legend = :topleft, margin=5Plots.mm, 
            legendfontsize=16, c=:red)
plot!([0], showaxis = false, grid = false, label = "Elastance of chamber walls", c=:gray30, ls=:dash);   
savefig(legend, "plots/chambers_volume_elastance_legend")


## PRESSURE, VOLUME AND ELASTANCE OF HEART CHAMBERS ##

# Left Ventricle
p_lv = plot(t, [sol[circ_sys.LV.p], sol[circ_sys.LV.V]], c=[:blue :red], xlabel="t (s)", ylabel="Pressure (mmHg) / Volume (ml)", 
            yguidefontsize=12, ytickfontsize=12, xguidefontsize=12, xtickfontsize=12, framestyle=:zerolines, 
            title = "Left Ventricle", margin=5Plots.mm, label = "")
plot!(twinx(), t, e_lv, c=:gray, ls=:dash, grid=true, xlabel="", ylabel="Elastance (mmHg/mL)", label="",  yguidefontsize=12, 
        ytickfontsize=12);

# Right Ventricle
p_rv = plot(t, [sol[circ_sys.RV.p], sol[circ_sys.RV.V]], c=[:blue :red], xlabel="t (s)", ylabel="Pressure (mmHg) / Volume (ml)", 
            yguidefontsize=12, ytickfontsize=12, xguidefontsize=12, xtickfontsize=12, framestyle=:zerolines, 
            title = "Right Ventricle", margin=5Plots.mm, label = "")
plot!(twinx(), t, e_rv, c=:gray, ls=:dash, grid=true, xlabel="", ylabel="Elastance (mmHg/mL)", label="",  yguidefontsize=12, 
        ytickfontsize=12);

# Left Atrium
p_la = plot(t, [sol[circ_sys.LA.p], sol[circ_sys.LA.V]], c=[:blue :red], xlabel="t (s)", ylabel="Pressure (mmHg) / Volume (ml)", 
            yguidefontsize=12, ytickfontsize=12, xguidefontsize=12, xtickfontsize=12, framestyle=:zerolines, 
            title = "Left Atrium", margin=5Plots.mm, label = "")
plot!(twinx(), t, e_la, c=:gray, ls=:dash, grid=true, xlabel="", ylabel="Elastance (mmHg/mL)", label="",  yguidefontsize=12, 
        ytickfontsize=12);

# Right Atrium
p_ra = plot(t, [sol[circ_sys.RA.p], sol[circ_sys.RA.V]], c=[:blue :red], xlabel="t (s)", ylabel="Pressure (mmHg) / Volume (ml)", 
            yguidefontsize=12, ytickfontsize=12, xguidefontsize=12, xtickfontsize=12, framestyle=:zerolines, 
            title = "Right Atrium", margin=5Plots.mm, label = "")
plot!(twinx(), t, e_ra, c=:gray, ls=:dash, grid=true, xlabel="", ylabel="Elastance (mmHg/mL)", label="",  yguidefontsize=12, 
        ytickfontsize=12);

# Adding four plots to one figure
fig = plot(p_lv, p_rv, p_la, p_ra, layout=(2,2), margin=5Plots.mm, legend = false)
plot!(size=(1200,700))
savefig(fig, "plots/chambers_pressure_volume_elastance")

# Create legend
legend = plot([0 0], showaxis = false, grid = false, label = ["Pressure in chamber" "Volume of blood in chamber"], 
        legend = :topleft,margin=5Plots.mm, legendfontsize=16, c=[:blue :red])
plot!([0], showaxis = false, grid = false, label = "Elastance of chamber walls", c=:gray30, ls=:dash);
savefig(legend, "plots/chambers_pressure_volume_elastance_legend")    


## PRESSURE ACROSS AND FLOW THROUGH HEART VALVES ##

# Mitral Valve
# creating y axis limits to ensure zero is at the same place on both y axes
a, b, c, d = extrema(sol[circ_sys.MV.Δp])..., extrema(sol[circ_sys.MV.q])...
α = -b/a
c * d ≥ 0 && (c = -d/α)
β = -d/c
α > β ? d = -α*c : c = -d/α

p_mv = plot(t, sol[circ_sys.MV.Δp], ylims=(a,b), c=:blue, label="Pressure across valve", xlabel="t (s)", ylabel="Pressure (mmHg)", 
            yguidefontsize=12, ytickfontsize=12, xguidefontsize=12, xtickfontsize=12, framestyle=:zerolines, title = "Mitral Valve", 
            margin=5Plots.mm, ls=:dash)
plot!((0, NaN), c=:red, label="Flow through valve")
plot!(twinx(), t, sol[circ_sys.MV.q], ylims=(c,d), c=:red, grid=true, xlabel="", ylabel="Flow rate (mL/s)", label="",  
        yguidefontsize=12, ytickfontsize=12)

# add tick marks on the left y axis
y = [25.0, 50.0, 75.0, 100.0];

for i in 1:length(y)
    plot!([0.0,0.01], [y[i],y[i]], c=:black, label="")
end

vline!([0.0], c=:black, label="")

# Tricuspid Valve
# creating y axis limits to ensure zero is at the same place on both y axes
a, b, c, d = extrema(sol[circ_sys.TV.Δp])..., extrema(sol[circ_sys.TV.q])...
α = -b/a
c * d ≥ 0 && (c = -d/α)
β = -d/c
α > β ? d = -α*c : c = -d/α

p_tv = plot(t, sol[circ_sys.TV.Δp], ylims=(a,b), c=:blue, label="Pressure across valve", xlabel="t (s)", ylabel="Pressure (mmHg)", 
            yguidefontsize=12, ytickfontsize=12, xguidefontsize=12, xtickfontsize=12, framestyle=:zerolines, title = "Tricuspid Valve", 
            margin=5Plots.mm, ls=:dash)
plot!((0, NaN), c=:red,label="Flow through valve")
plot!(twinx(), t, sol[circ_sys.TV.q], ylims=(c,d), c=:red, grid=true, xlabel="", ylabel="Flow rate (mL/s)", label="",  
        yguidefontsize=12, ytickfontsize=12)

# add tick marks on the left y axis
y = [10.0, 20.0, 30.0];

for i in 1:length(y)
    plot!([0.0,0.01], [y[i],y[i]], c=:black, label="")
end

vline!([0.0], c=:black, label="")

# Aortic Valve
# creating y axis limits to ensure zero is at the same place on both y axes
a, b, c, d = extrema(sol[circ_sys.AV.Δp])..., extrema(sol[circ_sys.AV.q])...
α = -b/a
c * d ≥ 0 && (c = -d/α)
β = -d/c
α > β ? d = -α*c : c = -d/α

p_av = plot(t, sol[circ_sys.AV.Δp], ylims=(a,b), c=:blue, label="Pressure across valve", xlabel="t (s)", ylabel="Pressure (mmHg)", 
            yguidefontsize=12, ytickfontsize=12, xguidefontsize=12, xtickfontsize=12, framestyle=:zerolines, title = "Aortic Valve", 
            margin=5Plots.mm, ls=:dash)
plot!((0, NaN), c=:red,label="Flow through valve")
plot!(twinx(), t, sol[circ_sys.AV.q], ylims=(c,d), c=:red, grid=true, xlabel="", ylabel="Flow rate (mL/s)", label="",  
        yguidefontsize=12, ytickfontsize=12)

# add tick marks on the left y axis
y = [25.0, 50.0, 75.0, 100.0];

for i in 1:length(y)
    plot!([0.0,0.01], [y[i],y[i]], c=:black, label="")
end

vline!([0.0], c=:black, label="")

# Pulmonary Valve
# creating y axis limits to ensure zero is at the same place on both y axes
a, b, c, d = extrema(sol[circ_sys.PV.Δp])..., extrema(sol[circ_sys.PV.q])...
α = -b/a
c * d ≥ 0 && (c = -d/α)
β = -d/c
α > β ? d = -α*c : c = -d/α

p_pv = plot(t, sol[circ_sys.PV.Δp], ylims=(a,b), c=:blue, label="Pressure across valve", xlabel="t (s)", ylabel="Pressure (mmHg)", 
            yguidefontsize=12, ytickfontsize=12, xguidefontsize=12, xtickfontsize=12, framestyle=:zerolines, title = "Pulmonary Valve", 
            margin=5Plots.mm, ls=:dash)
plot!((0, NaN), c=:red,label="Flow through valve")
plot!(twinx(), t, sol[circ_sys.PV.q], ylims=(c,d), c=:red, grid=true, xlabel="", ylabel="Flow rate (mL/s)", label="",  
        yguidefontsize=12, ytickfontsize=12)

# add tick marks on the left y axis
y = [10.0, 20.0, 30.0];

for i in 1:length(y)
    plot!([0.0,0.01], [y[i],y[i]], c=:black, label="")
end

vline!([0.0], c=:black, label="")

# Adding four plots to one figure
fig = plot(p_mv, p_tv, p_av, p_pv, layout=(2,2), margin=5Plots.mm, legend = false)
plot!(size=(1200,700))
savefig(fig, "plots/valves_pressure_flow")

# Create legend
legend = plot([0], showaxis = false, grid = false, label = "Pressure across valve" , legend = :topleft, margin=5Plots.mm, 
            legendfontsize=16, c=:blue, ls=:dash)
plot!([0], showaxis = false, grid = false, label = "Flow through valve", c=:red);
savefig(legend, "plots/valves_pressure_flow_legend")   