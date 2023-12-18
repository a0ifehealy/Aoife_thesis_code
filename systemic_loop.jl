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

# OnePort
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

##########################################
# BLOOD FLOW FUNCTION (SIMPLE SINE WAVE) #
##########################################

A = 314     # mL/s
B = 1;      # s

function flow(t)
    bloodflow =  A*sin(2*π*t/B)
    return max(0,bloodflow)
end

#################################################################
# DEFINING PARAMETERS THAT WILL BE USED IN SYSTEMIC CIRCULATION #
#################################################################

L_ao = 6.2e-5; C_ao = 0.08; R_ao = 0.003;
L_ar = 0.0017; C_ar = 1.6; R_ar = 0.05;
R_at = 0.5;
R_cp = 0.52;
C_vn = 20.5; R_vn = 0.075;

##################
# BUILDING MODEL #
##################

@mtkmodel systemic_circulation_model begin
    @components begin
        source = DrivenCurrent(Q = 1.0)
        ground = Ground()
        # Aorta
        resistor_ao = Resistor(R = R_ao)
        capacitor_ao = Capacitor(C = C_ao)
        inductor_ao = Inductor(L = L_ao)
        # Arteries
        resistor_ar = Resistor(R = R_ar)
        capacitor_ar = Capacitor(C = C_ar)
        inductor_ar = Inductor(L = L_ar)
        # Arterioles
        resistor_at = Resistor(R = R_at)
        # Capilliaries
        resistor_cp = Resistor(R = R_cp)
        # Veins
        resistor_vn = Resistor(R = R_vn)
        capacitor_vn = Capacitor(C = C_vn)
    end
    @equations begin
        connect(source.out, inductor_ao.in, capacitor_ao.in)
        connect(capacitor_ao.out, ground.g)
        connect(inductor_ao.out, resistor_ao.in)
        connect(resistor_ao.out, inductor_ar.in, capacitor_ar.in)
        connect(capacitor_ar.out, ground.g)
        connect(inductor_ar.out, resistor_ar.in)
        connect(resistor_ar.out, resistor_at.in)
        connect(resistor_at.out, resistor_cp.in)
        connect(resistor_cp.out, resistor_vn.in, capacitor_vn.in)
        connect(resistor_vn.out, capacitor_vn.out, ground.g)
        connect(source.in, ground.g)
    end
end

@mtkbuild systemic_circulation = systemic_circulation_model()

###############
# SOLVING ODE #
###############

u0 = [
    systemic_circulation.capacitor_ao.Δp => -78,
    systemic_circulation.inductor_ao.q => 50,
    systemic_circulation.capacitor_ar.Δp => -50,
    systemic_circulation.inductor_ar.q => 50,
    systemic_circulation.capacitor_vn.Δp => -20
]

prob = ODEProblem(systemic_circulation, u0, (0, 25.0))
sol= solve(prob, RK4(), reltol=1e-6)
plot(sol, idxs = systemic_circulation.source.Δp, tspan = (20, 25))