using ModelingToolkit
using DifferentialEquations

export Pin, Ground, OnePort, Resistor, Capacitor, Compliance, Inductance, DrivenFlow, HeartChamber, HeartValve, CRL, CR

@variables t
D = Differential(t)

# Pin

@connector function Pin(; name)
    sts = @variables p(t) = 1.0 q(t) = 1.0 [connect = Flow]
    ODESystem(Equation[], t, sts, []; name=name)
end

# Ground

@component function Ground(; name, P=0.0)
    @named g = Pin()
    ps = @parameters P = P
    eqs = [g.p ~ P]
    compose(ODESystem(eqs, t, [], ps; name=name), g)
end

# OnePort

@component function OnePort(; name)
    @named in = Pin()
    @named out = Pin()
    sts = @variables Δp(t) = 0.0 q(t) = 0.0
    eqs = [
            Δp ~ out.p - in.p
            0 ~ in.q + out.q
            q ~ in.q
    ]
    compose(ODESystem(eqs, t, sts, []; name=name), in, out)
end

# Resistor

@component function Resistor(; name, R=1.0)
    @named oneport = OnePort()
    @unpack Δp, q = oneport
    ps = @parameters R = R
    eqs = [
            Δp ~ -q * R
    ]
    extend(ODESystem(eqs, t, [], ps; name=name), oneport)
end

# Capacitor

@component function Capacitor(; name, C=1.0)
    @named oneport = OnePort()
    @unpack Δp, q = oneport
    ps = @parameters C = C
    eqs = [
            D(Δp) ~ -q / C
    ]
    extend(ODESystem(eqs, t, [], ps; name=name), oneport)
end

# Compliance

@component function Compliance(; name, V₀=0.0, C=1.0, p₀=0.0)
    @named in = Pin()
    @named out = Pin()

    sts = @variables begin
            V(t) = V₀
            p(t) = 0.0
    end

    ps = @parameters begin
            V₀ = V₀
            C = C
            p_rel = p₀
    end

    D = Differential(t)

    eqs = [
            0 ~ in.p - out.p
            p ~ in.p
            p ~ (V - V₀) / C + p_rel
            D(V) ~ in.q + out.q
    ]

    compose(ODESystem(eqs, t, sts, ps; name=name), in, out)
end

# Inductor

@component function Inductance(; name, L=1.0)
    @named oneport = OnePort()
    @unpack Δp, q = oneport
    ps = @parameters L = L
    D = Differential(t)
    eqs = [
            D(q) ~ -Δp / L
    ]
    extend(ODESystem(eqs, t, [], ps; name=name), oneport)
end

# Driven Current

@component function DrivenFlow(; name, Q=1.0)
    @named oneport = OnePort()
    @unpack q = oneport
    ps = @parameters Q = Q
    eqs = [
            q ~ Q * flow(t)
    ]
    extend(ODESystem(eqs, t, [], ps; name=name), oneport)
end

# Heart Chamber
@component function HeartChamber(; name, V₀, p₀, E_min, E_max, T, T_es, T_ep, Eshift=0.0)
    @named in = Pin()
    @named out = Pin()
    sts = @variables V(t) = 0.0 p(t) = 0.0
    ps = @parameters V₀ = V₀ p₀ = p₀ E_min = E_min E_max = E_max T = T T_es = T_es T_ep = T_ep Eshift = Eshift

    D = Differential(t)
    E = Elastance(t, E_min, E_max, T, T_es, T_ep, Eshift)

    p_rel = p₀

    eqs = [
            0 ~ in.p - out.p
            p ~ in.p
            p ~ (V - V₀) * E + p_rel
            D(V) ~ in.q + out.q
        ]

    compose(ODESystem(eqs, t, sts, ps; name=name), in, out)
end

# Valve

@component function HeartValve(; name, CQ=1.0)
    @named oneport = OnePort()
    @unpack Δp, q = oneport
    ps = @parameters CQ = CQ
    eqs = [
            q ~ (Δp < 0) * CQ * sqrt(abs(Δp))
    ]
    extend(ODESystem(eqs, t, [], ps; name=name), oneport)
end

# CRL

@component function CRL(; name, C=1.0, R=1.0, L=1.0)
    @named in = Pin()
    @named out = Pin()

    sts = @variables Δp(t) = 0.0 q(t) = 0.0
    ps = []

    @named C = Compliance(C=C)
    @named R = Resistor(R=R)
    @named L = Inductance(L=L)

    eqs = [
            Δp ~ out.p - in.p
            q ~ in.q
            connect(in, C.in)
            connect(C.out, R.in)
            connect(R.out, L.in)
            connect(L.out, out)
    ]

    compose(ODESystem(eqs, t, sts, ps; name=name), in, out, C, R, L)
end

# CR

@component function CR(; name, R=1.0, C=1.0)
    @named in = Pin()
    @named out = Pin()

    sts = @variables Δp(t) = 0.0 q(t) = 0.0
    ps = []

    @named R = Resistor(R=R)
    @named C = Compliance(C=C)

    eqs = [
            Δp ~ out.p - in.p
            q ~ in.q
            connect(in, C.in)
            connect(C.out, R.in)
            connect(R.out, out)
    ]

    compose(ODESystem(eqs, t, sts, ps; name=name), in, out, R, C)
end