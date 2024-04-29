using ModelingToolkit
using DifferentialEquations

export  Pin, Ground, OnePort, Resistor, Capacitor, Compliance, Inductance, 
        DrivenFlow, HeartChamber, HeartValve, CRL, CR

@variables t
D = Differential(t)

# Pin
"""
Pin(; name)

Defines a node, at which pressure and flow rate is continuous
"""
@connector function Pin(; name)
    sts = @variables p(t) = 1.0 q(t) = 1.0 [connect = Flow]
    ODESystem(Equation[], t, sts, []; name=name)
end

# Ground
"""
Ground(; name)

Defines the ground node, or a reference node. Pressure defaults to zero.
"""
@component function Ground(; name, P=0.0)
    @named g = Pin()
    ps = @parameters P = P
    eqs = [g.p ~ P]
    compose(ODESystem(eqs, t, [], ps; name=name), g)
end

# OnePort
"""
OnePort(; name)

Defines a generic 'electrical' element by applying Kirchhoff's circuit 
laws.

Parameters calculated:
del_p      Pressure drop across element (mmHg)
q       Blood flow through element (ml/s)
"""
@component function OnePort(; name)
    @named in = Pin()
    @named out = Pin()
    sts = @variables del_p(t) = 0.0 q(t) = 0.0
    eqs = [
            del_p ~ out.p - in.p
            0 ~ in.q + out.q
            q ~ in.q
    ]
    compose(ODESystem(eqs, t, sts, []; name=name), in, out)
end

# Resistor
"""
Resistor(; name, R=1.0)

Defines a resistor by applying Ohm's law to a generic element.

Arguments:
R       Resistance of blood vessel (mmHg*s/ml)

Parameters calculated:
del_p      Pressure drop across resistance element (mmHg)
q       Blood flow through resistance element (ml/s)
"""
@component function Resistor(; name, R=1.0)
    @named oneport = OnePort()
    @unpack del_p, q = oneport
    ps = @parameters R = R
    eqs = [
            del_p ~ -q * R
    ]
    extend(ODESystem(eqs, t, [], ps; name=name), oneport)
end

# Capacitor
"""
Capacitor(; name, C=1.0)

Defines a capacitor by applying the capacitor's pressure-flow 
relationship to a generic element.

Arguments:
C       Compliance of blood vessel (ml/mmHg)

Parameters calculated:
del_p      Pressure drop across capacitance element (mmHg)
q       Blood flow through capacitance element (ml/s)
"""
@component function Capacitor(; name, C=1.0)
    @named oneport = OnePort()
    @unpack del_p, q = oneport
    ps = @parameters C = C
    eqs = [
            D(del_p) ~ -q / C
    ]
    extend(ODESystem(eqs, t, [], ps; name=name), oneport)
end

# Compliance
"""
Compliance(; name, V_0=0.0, C=1.0, p_0=0.0)

Defines a compliant vessel.

Arguments:
V_0      Stress-free (zero pressure) volume (ml)
C       Compliance of blood vessel (ml/mmHg)
p_0      Offset pressure value (mmHg)

Parameters calculated:
p       Pressure in compliance vessel (mmHg)
V       Volume in compliance vessel (ml)
"""
@component function Compliance(; name, V_0=0.0, C=1.0, p_0=0.0)
    @named in = Pin()
    @named out = Pin()

    sts = @variables begin
            V(t) = V_0
            p(t) = 0.0
    end

    ps = @parameters begin
            V_0 = V_0
            C = C
            p_rel = p_0
    end

    eqs = [
            0 ~ in.p - out.p
            p ~ in.p
            p ~ ((V - V_0) / C) + p_rel
            D(V) ~ in.q + out.q
    ]

    compose(ODESystem(eqs, t, sts, ps; name=name), in, out)
end

# Inductor
"""
Inductance(; name, L=1.0)

Defines a capacitor by applying the inductors's pressure-flow 
relationship to a generic element to represent inertia.

Arguments:
L       Inductance of blood vessel (mmHg*s^2/ml)

Parameters calculated:
del_p      Pressure drop across capacitance element (mmHg)
q       Blood flow through capacitance element (ml/s)
"""
@component function Inductance(; name, L=1.0)
    @named oneport = OnePort()
    @unpack del_p, q = oneport
    ps = @parameters L = L
    D = Differential(t)
    eqs = [
            D(q) ~ -del_p / L
    ]
    extend(ODESystem(eqs, t, [], ps; name=name), oneport)
end

# Driven Current
"""
DrivenFlow(; name, Q=1.0)

Defines a blood flow source.

Arguments:
Q       Magnitude of blood flow  (ml/s)

Parameters calculated:
q       Blood flow through source element (ml/s)
"""
@component function DrivenFlow(; name, Q=1.0)
    @named oneport = OnePort()
    @unpack q = oneport
    ps = @parameters Q = Q
    eqs = [
            q ~ Q * flow(t)
    ]
    extend(ODESystem(eqs, t, [], ps; name=name), oneport)
end

"""
HeartChamber(; name, V_0, p_0, E_min, E_max, T, T_es, T_ep, Eshift=0.0)

Defines a chamber of the heart

Arguments:
V_0      Stress-free (zero pressure) volume (ml)
p_0      Pressure offset value (mmHg)
E_min   Minimum elastance (mmHg/ml)
E_max   Maximum elastance (mmHg/ml)
T       Period of cardiac cycle (s)
T_es    End systolic time (s)
T_ep    Peak systolic time (s)
Eshift  Time shift of contraction (s)

Parameters calculated:
p       Pressure in chamber (mmHg)
V       Volume in chamber (ml)
"""
# Heart Chamber
@component function HeartChamber(; name, V_0, p_0, E_min, E_max, T, T_es, 
    T_ep, Eshift=0.0)

    @named in = Pin()
    @named out = Pin()
    sts = @variables V(t) = 0.0 p(t) = 0.0
    ps = @parameters (V_0 = V_0, p_0 = p_0, E_min = E_min, E_max = E_max, 
                    T = T, T_es = T_es, T_ep = T_ep, Eshift = Eshift)

    D = Differential(t)
    E = Elastance(t, E_min, E_max, T, T_es, T_ep, Eshift)

    p_rel = p_0

    eqs = [
            0 ~ in.p - out.p
            p ~ in.p
            p ~ ((V - V_0) * E) + p_rel
            D(V) ~ in.q + out.q
        ]

    compose(ODESystem(eqs, t, sts, ps; name=name), in, out)
end

# Valve
"""
HeartValve(; name, CQ=1.0)

Defines the heart valves as orifice valves

Arguments:
CQ      Flow coefficient (ml/(s*mmHg^0.5))

Parameters calculated:
del_p      Pressure drop across capacitance element (mmHg)
q       Blood flow through capacitance element (ml/s)
"""
@component function HeartValve(; name, CQ=1.0)
    @named oneport = OnePort()
    @unpack del_p, q = oneport
    ps = @parameters CQ = CQ
    eqs = [
            q ~ (del_p < 0) * CQ * sqrt(abs(del_p))
    ]
    extend(ODESystem(eqs, t, [], ps; name=name), oneport)
end

# CRL
"""
CRL(; name, C=1.0, R=1.0, L=1.0)

Defines compliace, resistor, inductance subsystem

Arguments:
C       Blood vessel compliace (ml/mmHg)
R       Blood vessel resistance (mmHg*s/ml)
L       Blood vessel inductance (mmHg*s^2/ml)

Parameters calculated:
del_p      Pressure drop across subsystem (mmHg)
q       Blood flow through subsystem element (ml/s)
"""
@component function CRL(; name, C=1.0, R=1.0, L=1.0)
    @named in = Pin()
    @named out = Pin()

    sts = @variables del_p(t) = 0.0 q(t) = 0.0
    ps = []

    @named C = Compliance(C=C)
    @named R = Resistor(R=R)
    @named L = Inductance(L=L)

    eqs = [
            del_p ~ out.p - in.p
            q ~ in.q
            connect(in, C.in)
            connect(C.out, R.in)
            connect(R.out, L.in)
            connect(L.out, out)
    ]

    compose(ODESystem(eqs, t, sts, ps; name=name), in, out, C, R, L)
end

# CR
"""
CR(; name, C=1.0, R=1.0)

Defines compliace, resistor subsystem

Arguments:
C       Blood vessel compliace (ml/mmHg)
R       Blood vessel resistance (mmHg*s/ml)

Parameters calculated:
del_p      Pressure drop across subsystem (mmHg)
q       Blood flow through subsystem element (ml/s)
"""
@component function CR(; name, R=1.0, C=1.0)
    @named in = Pin()
    @named out = Pin()

    sts = @variables del_p(t) = 0.0 q(t) = 0.0
    ps = []

    @named R = Resistor(R=R)
    @named C = Compliance(C=C)

    eqs = [
            del_p ~ out.p - in.p
            q ~ in.q
            connect(in, C.in)
            connect(C.out, R.in)
            connect(R.out, out)
    ]

    compose(ODESystem(eqs, t, sts, ps; name=name), in, out, R, C)
end