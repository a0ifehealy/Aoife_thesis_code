using DifferentialEquations
using Plots

# Define parameters that will be used in Windkessel models
A = 100         # Blood flow rate mL/s
B = 60/72       # Period of cardiac cycle (s)
C = 1.1         # Compliance (ml/mmHg)
R = 0.9         # Distal resistance (mmHg*s/ml)
r = 0.05        # Proximal Resistance (mmHg*s/ml)
L = 0.008       # Blood inertia (mmHg*s^2/ml)
tspan = (0,5)   # Time period over which ODEs are solved
P0 = 78         # Pressure (in aorta?) when t = 0

# Blood flow function (simple sine wave)
function flow(t)
    return A*sin(2*π*t/B)+A
end

function flow_dot(t)
    return (2*A*π/B)*cos(2*π*t/B)
end

function flow_dotdot(t)
    return -(4*A*(π^2)/(B^2))*sin(2*π*t/B)
end

# 2 element Windkessel
function pressure_2_element(P,p,t)
    return flow(t)/C - P/(C*R)
end

prob_2_elem = ODEProblem(pressure_2_element,P0,tspan)
sol_2_elem = solve(prob_2_elem,RK4(), reltol=1e-6);

# 3 element Windkessel
function pressure_3_element(P,p,t)
    I = flow(t)
    Idot = flow_dot(t)
    return (I/C)*(1+r/R) + r*Idot - P/(R*C)
end

prob_3_elem = ODEProblem(pressure_3_element,P0,tspan)
sol_3_elem = solve(prob_3_elem,RK4(), reltol=1e-6);

# 4 element Windkessel
function pressure_4_element(P,p,t)
    I = flow(t)
    Idot = flow_dot(t)
    Idotdot = flow_dotdot(t)
    return L*Idotdot + (L/(C*R)+r)*Idot + (1/C+r/(R*C))*I - P/(C*R)
end

prob_4_elem = ODEProblem(pressure_4_element,P0,tspan)
sol_4_elem = solve(prob_4_elem,RK4(), reltol=1e-6);

# Plotting
plot(sol_2_elem, labels="2 element", ylabel="Voltage", xlabel="Time")
plot!(sol_3_elem, labels="3 element")
plot!(sol_4_elem, labels="4 element")