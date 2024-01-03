using ModelingToolkit
using DifferentialEquations
using Plots

@variables t
D = Differential(t)
include("model_parameters.jl")
include("circulatory_system_components.jl")

######################
# BLOODFLOW FUNCTION #
######################

# Define bloodflow into aorta as a simple sine wave
function flow(t)
    bloodflow = A_sys_pul*sin(2*π*t/T)
    return max(0,bloodflow)
end

# Store bloodflow in aorta data for 3 cycles
t1 = 0:0.01:3*T
bloodflow = []
for t in t1
    append!(bloodflow,flow(t))
end

################################################################
# DEFINE COMPONENTS USED IN BOTH SYSTEMIC AND PULMONAARY LOOPS #
################################################################

# Source of blood flow
@named source = DrivenFlow(Q = 1.0)

# Referenece pressure
@named ground = Ground()

#############################################
# DEFINE COMPONENTS OF SYSTEMIC CIRCULATION #
#############################################

# Systemic Aortic Sinus
@named R_sas = Resistor(R = Rsas)
@named C_sas = Capacitor(C = Csas)
@named L_sas = Inductance(L = Lsas)

# Systemic Arteries
@named R_sat = Resistor(R = Rsat)
@named C_sat = Capacitor(C = Csat)
@named L_sat = Inductance(L = Lsat)

# Systemic Arterioles
@named R_sar = Resistor(R = Rsar)

# Systemic Capillaries
@named R_scp = Resistor(R = Rscp)

# Systemic Arteries
@named R_svn = Resistor(R = Rsvn)
@named C_svn = Capacitor(C = Csvn)

###################################
# COMPOSE AND SOLVE SYSTEMIC LOOP #
###################################

# Connect components to form circuit
sys_eqs = [
    connect(source.out, L_sas.in, C_sas.in)
    connect(C_sas.out, ground.g)
    connect(L_sas.out, R_sas.in)
    connect(R_sas.out, L_sat.in, C_sat.in)
    connect(C_sat.out, ground.g)
    connect(L_sat.out, R_sat.in)
    connect(R_sat.out, R_sar.in)
    connect(R_sar.out, R_scp.in)
    connect(R_scp.out, R_svn.in, C_svn.in)
    connect(R_svn.out, C_svn.out, ground.g)
    connect(source.in, ground.g)
];

# Create ODE system
@named _sys_model = ODESystem(sys_eqs, t)
@named sys_model = compose(_sys_model, [source, L_sas, C_sas, R_sas, L_sat, C_sat, R_sat, R_sar, R_scp, R_svn, C_svn, ground])
sys_circ = structural_simplify(sys_model)

#Initial conditions for ODEs
u0_sys = [
    sys_circ.C_sas.Δp => -pt0sas,
    sys_circ.L_sas.q => qt0sas,
    sys_circ.C_sat.Δp => -pt0sat,
    sys_circ.L_sat.q => qt0sat,
    sys_circ.C_svn.Δp => -pt0svn
]

# Solve system of ODEs
prob = ODEProblem(sys_circ, u0_sys, sys_pul_tspan)
sol_sys = solve(prob, RK4(), reltol=1e-6)

#############################################
# DEFINE COMPONENTS OF PULMONARY CIRCULATION #
#############################################

# Systemic Aortic Sinus
@named R_pas = Resistor(R = Rsas)
@named C_pas = Capacitor(C = Cpas)
@named L_pas = Inductance(L = Lpas)

# Systemic Arteries
@named R_pat = Resistor(R = Rpat)
@named C_pat = Capacitor(C = Cpat)
@named L_pat = Inductance(L = Lpat)

# Systemic Arterioles
@named R_par = Resistor(R = Rpar)

# Systemic Capillaries
@named R_pcp = Resistor(R = Rpcp)

# Systemic Arteries
@named R_pvn = Resistor(R = Rpvn)
@named C_pvn = Capacitor(C = Cpvn)

###################################
# COMPOSE AND SOLVE PULMONARY LOOP #
###################################

# Connect components to form circuit
pul_eqs = [
    connect(source.out, L_pas.in, C_pas.in)
    connect(C_pas.out, ground.g)
    connect(L_pas.out, R_pas.in)
    connect(R_pas.out, L_pat.in, C_pat.in)
    connect(C_pat.out, ground.g)
    connect(L_pat.out, R_pat.in)
    connect(R_pat.out, R_par.in)
    connect(R_par.out, R_pcp.in)
    connect(R_pcp.out, R_pvn.in, C_pvn.in)
    connect(R_pvn.out, C_pvn.out, ground.g)
    connect(source.in, ground.g)
];

# Create ODE system
@named _pul_model = ODESystem(pul_eqs, t)
@named pul_model = compose(_pul_model, [source, L_pas, C_pas, R_pas, L_pat, C_pat, R_pat, R_par, R_pcp, R_pvn, C_pvn, ground])
pul_circ = structural_simplify(pul_model)

#Initial conditions for ODEs
u0_pul = [
    pul_circ.C_pas.Δp => -pt0pas,
    pul_circ.L_pas.q => qt0pas,
    pul_circ.C_pat.Δp => -pt0pat,
    pul_circ.L_pat.q => qt0pat,
    pul_circ.C_pvn.Δp => -pt0pvn
]

# Solve system of ODEs
prob = ODEProblem(pul_circ, u0_pul, sys_pul_tspan)
sol_pul = solve(prob, RK4(), reltol=1e-6)

####################
# PLOTTING RESULTS #
####################

# Re-indexing results so that time starts after the solution has stabilised
function find_n(sol, t_min)
    n = 1
    for t in sol.t
        if t>t_min
            break
        end
        n+=1
    end
    return n
end

n_sys = find_n(sol_sys, sys_pul_reindex_t);
t_sys = sol_sys.t[n_sys:end] .- sol_sys.t[n_sys]
y_sys = sol_sys[source.Δp][n_sys:end];

n_pul = find_n(sol_pul, sys_pul_reindex_t);
t_pul = sol_pul.t[n_pul:end] .- sol_pul.t[n_pul]
y_pul = sol_pul[source.Δp][n_pul:end];

# Plot pressure in aortic sinuses for systemic and pulmonary circulation
p1 = plot(t_sys, y_sys, ylabel="Aortic Sinus Pressure (mmHg)", xlabel="t (s)",legend = false, margin=5Plots.mm, 
            yguidefontsize=20, ytickfontsize=16, xguidefontsize=20, xtickfontsize=16, lw=2, xlims=(0,3.0))
plot!(t_pul, y_pul, lw=2)
plot!(twinx(), t1, bloodflow, c=:gray, ls=:dash, xlims=(0,3.0), label="", ylabel="Flow rate (ml/s)",
        yguidefontsize=20, ytickfontsize=16, lw = 2)
plot!(size=(1200,700))
savefig(p1, "plots/systemic_pulmonary")

# Create legend (needs to be done seperately to avoid bug when plotting using twinx())
legend = plot([0 0], showaxis = false, grid = false, label = ["Systemic Aortic Sinus Pressure" "Pulmonary Aortic Sinus Pressure"], 
        legend = :topleft, margin=5Plots.mm, legendfontsize=16,)
plot!([0], showaxis = false, grid = false, label = "Blood flow into Aortic Sinuses", c=:gray, ls=:dash) 
savefig(legend, "plots/systemic_pulmonary_legend")       
