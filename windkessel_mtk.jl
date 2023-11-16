using ModelingToolkit
using Plots
using DifferentialEquations

@variables t
D = Differential(t)

#############################
# DEFINING BASIC COMPONENTS #
#############################

# Pin (node)
@connector Pin begin
    p(t)
    q(t), [connect = Flow]
end

# Ground
@mtkmodel Ground begin
    @components begin
        g = Pin()
    end
    @equations begin
        g.p ~ 0
    end
end

# OnePort (generic element)
@mtkmodel OnePort begin
    @components begin
        in = Pin()
        out = Pin()
    end
    @variables begin
        Δp(t)
        q(t)
    end
    @equations begin
        Δp ~ out.p - in.p
        0 ~ in.q + out.q
        q ~ in.q
    end
end

# Resistor
@mtkmodel Resistor begin
    @extend OnePort()
    @parameters begin
        R = 1.0 # Sets the default resistance
    end
    @equations begin
        Δp ~ -q * R
    end
end

# Capacitor
@mtkmodel Capacitor begin
    @extend OnePort()
    @parameters begin
        C = 1.0
    end
    @equations begin
        D(Δp) ~ -q / C
    end
end

# Inductor
@mtkmodel Inductor begin
    @extend OnePort()
    @parameters begin
        L = 1.0
    end
    @equations begin
        D(q) ~ -Δp / L
    end
end

# Driven Current
@mtkmodel DrivenCurrent begin
    @extend OnePort()
    @parameters begin
        Q = 1.0
    end
    @equations begin
        q ~ Q * flow(t)
    end
end

##############################################################
# DEFINING PARAMETERS THAT WILL BE USED IN WINDKESSEL MODELS #
##############################################################

A = 100         # Blood flow rate mL/s
B = 60/72       # Period of cardiac cycle (s)
C = 1.1         # Compliance (ml/mmHg)
R = 0.9         # Distal resistance (mmHg*s/ml)
r = 0.05        # Proximal Resistance (mmHg*s/ml)
L = 0.008       # Blood inertia (mmHg*s^2/ml)
tspan = (0,5)   # Time period over which ODEs are solved
P0 = -78         # Pressure (in aorta?) when t = 0

##########################################
# BLOOD FLOW FUNCTION (SIMPLE SINE WAVE) #
##########################################

function flow(t)
    return A*sin(2*π*t/B)+A
end

########################
# 2 ELEMENT WINDKESSEL #
########################

@mtkmodel WK2model begin
    @components begin
        resistor = Resistor(R = R)
        capacitor = Capacitor(C = C)
        source = DrivenCurrent(Q = 1.0)
        ground = Ground()
    end
    @equations begin
        connect(source.out, resistor.in, capacitor.in)
        connect(resistor.out, source.in, capacitor.out, ground.g)
    end
end

@mtkbuild wk2model = WK2model()
u0 = [
    wk2model.capacitor.Δp => P0
]
prob = ODEProblem(wk2model, u0, tspan)
sol_2_elem = solve(prob, RK4(), reltol=1e-6)

########################
# 3 ELEMENT WINDKESSEL #
########################

@mtkmodel WK3model begin
    @components begin
        resistor_p = Resistor(R = R)
        resistor_c = Resistor(R = r)
        capacitor = Capacitor(C = C)
        source = DrivenCurrent(Q = 1.0)
        ground = Ground()
    end
    @equations begin
        connect(source.out, resistor_c.in)
        connect(resistor_c.out, resistor_p.in, capacitor.in)
        connect(resistor_p.out, capacitor.out, source.in, ground.g)
    end
end

@mtkbuild wk3model = WK3model()
u0 = [
    wk3model.capacitor.Δp => P0
]
prob = ODEProblem(wk3model, u0, tspan)
sol_3_elem = solve(prob, RK4(), reltol=1e-6)

########################
# 4 ELEMENT WINDKESSEL #
########################

@mtkmodel WK4model begin
    @components begin
        resistor_p = Resistor(R = R)
        resistor_c = Resistor(R = r)
        capacitor = Capacitor(C = C)
        inductor = Inductor(L=L)
        source = DrivenCurrent(Q = 1.0)
        ground = Ground()
    end
    @equations begin
        connect(source.out, resistor_c.in)
        connect(resistor_c.out, inductor.in)
        connect(inductor.out, resistor_p.in, capacitor.in)
        connect(resistor_p.out, capacitor.out, source.in, ground.g)
    end
end

@mtkbuild wk4model = WK4model()
u0 = [
    wk4model.capacitor.Δp => P0
]
prob = ODEProblem(wk4model, u0, tspan)
sol_4_elem = solve(prob, RK4(), reltol=1e-6)

############
# PLOTTING #
############

plot(sol_2_elem, idxs = wk2model.source.Δp, labels="2 element", xlabel="Time", ylabel="Pressure")
plot!(sol_3_elem, idxs = wk3model.source.Δp, labels="3 element")
plot!(sol_4_elem, idxs = wk4model.source.Δp, labels="4 element")