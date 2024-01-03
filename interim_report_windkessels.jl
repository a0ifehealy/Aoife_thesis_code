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
    bloodflow = wk_A_sine*sin(2*π*t/T-π/2) + wk_A_sine
    return bloodflow
end

# Store bloodflow in aorta data for 3 cycles
t1 = 0:0.01:3*T
bloodflow = []
for t in t1
    append!(bloodflow,flow(t))
end

############################################
# DEFINE COMPONENTS OF WINDKESSEL CIRCUITS #
############################################

# Peripheral resistance
@named resistor_p = Resistor(R = wk_rp)
# Characteristic impedance
@named resistor_c = Resistor(R = wk_rc)
# Compliance (capacitance)
@named capacitor = Capacitor(C = wk_c)
# Inertia (inductance)
@named inductor = Inductance(L = wk_l)
# Source of blood flow
@named source = DrivenFlow(Q = 1.0)
# Referenece pressure
@named ground = Ground()

########################
# 2 ELEMENT WINDKESSEL #
########################

# Connect components to form circuit
wk_2_eqs = [
    connect(source.out, resistor_p.in, capacitor.in)
    connect(resistor_p.out, source.in, capacitor.out, ground.g)
]

# Create ODE system
@named _wk2_model = ODESystem(wk_2_eqs, t)
@named wk2_model = compose(_wk2_model, [source, resistor_p, capacitor, ground])
wk2_sys = structural_simplify(wk2_model)

#Initial conditions for ODEs
u0 = [
    wk2_sys.capacitor.Δp => wk_capacitor_delp0
]

# Solve system of ODEs
prob = ODEProblem(wk2_sys, u0, wk_tspan)
sol_2_elem = solve(prob, RK4(), reltol=1e-6)

########################
# 3 ELEMENT WINDKESSEL #
########################

# Connect components to form circuit
wk_3_eqs = [
    connect(source.out, resistor_c.in)
    connect(resistor_c.out, resistor_p.in, capacitor.in)
    connect(resistor_p.out, capacitor.out, source.in, ground.g)
]

# Create ODE system
@named _wk3_model = ODESystem(wk_3_eqs, t)
@named wk3_model = compose(_wk3_model, [source, resistor_c, resistor_p, capacitor, ground])
wk3_sys = structural_simplify(wk3_model)

#Initial conditions for ODEs
u0 = [
    wk3_sys.capacitor.Δp => wk_capacitor_delp0
]

# Solve system of ODEs
prob = ODEProblem(wk3_sys, u0, wk_tspan)
sol_3_elem = solve(prob, RK4(), reltol=1e-6)

########################
# 4 ELEMENT WINDKESSEL #
########################

# Connect components to form circuit
wk_4_eqs = [
    connect(source.out, resistor_c.in)
    connect(resistor_c.out, inductor.in)
    connect(inductor.out, resistor_p.in, capacitor.in)
    connect(resistor_p.out, capacitor.out, source.in, ground.g)
]

# Create ODE system
@named _wk4_model = ODESystem(wk_4_eqs, t)
@named wk4_model = compose(_wk4_model, [source, resistor_c, resistor_p, capacitor, inductor, ground])
wk4_sys = structural_simplify(wk4_model)

#Initial conditions for ODEs
u0 = [
    wk4_sys.capacitor.Δp => wk_capacitor_delp0
]

# Solve system of ODEs
prob = ODEProblem(wk4_sys, u0, wk_tspan)
sol_4_elem = solve(prob, RK4(), reltol=1e-6)

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

n_2wk = find_n(sol_2_elem, wk_reindex_t);
t_sol2 = sol_2_elem.t[n_2wk:end] .- sol_2_elem.t[n_2wk]
y_sol2 = sol_2_elem[source.Δp][n_2wk:end];

n_3wk = find_n(sol_3_elem, wk_reindex_t);
t_sol3 = sol_3_elem.t[n_3wk:end] .- sol_3_elem.t[n_3wk]
y_sol3 = sol_3_elem[source.Δp][n_3wk:end];

n_4wk = find_n(sol_4_elem, wk_reindex_t);
t_sol4 = sol_4_elem.t[n_4wk:end] .- sol_4_elem.t[n_4wk]
y_sol4 = sol_4_elem[source.Δp][n_4wk:end];

# Plot the pressure in aorta from the 3 different windkessel models
windkessel_plot = plot(t_sol2, y_sol2, labels="2 Element WK", xlabel="t (s)", ylabel="Aortic Pressure (mmHg)", legend=false,
        margin=5Plots.mm, yguidefontsize=20, ytickfontsize=16, xguidefontsize=20, xtickfontsize=16, lw=2, xlims=(0.0,3.0))
plot!(t_sol3, y_sol3, labels="3 Element WK", lw = 2, xlims=(0.0,3.0))
plot!(t_sol4, y_sol4, labels="4 Element WK", lw = 2, xlims=(0.0,3.0))
plot!(twinx(), t1, bloodflow, c=:gray, ls=:dash, xlims=(0.0,3.0), label="", ylabel="Flow rate (ml/s)",
        yguidefontsize=20, ytickfontsize=16, lw = 2)
plot!((0, NaN), c=:gray, ls=:dash, label="Blood Flow (ml/s)")
plot!(size=(1200,700))
savefig(windkessel_plot, "plots/windkessel_simple_sine_input")

# Create legend (needs to be done seperately to avoid bug when plotting using twinx())
windkessel_legend = plot([0 0 0], showaxis = false, grid = false, 
            label = ["Aortic Pressure, 2 Element WK" "Aortic Pressure, 3 Element WK" "Aortic Pressure, 4 Element WK"], 
            legend = :topleft, margin=5Plots.mm, lw = 2, legendfontsize=16,)
plot!([0], showaxis = false, grid = false, label = "Blood flow into Aorta", c=:gray, ls=:dash)
savefig(windkessel_legend, "plots/windkessel_legend")        