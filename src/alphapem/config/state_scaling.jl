# -*- coding: utf-8 -*-

"""
Module StateScalingModule

Centralized scaling references for solver state variables.

Why scaling is introduced
-------------------------
The ODE state vector mixes variables with very different physical magnitudes:
gas concentrations are typically O(10), temperatures are O(10²), cathode
overpotential is O(10⁻¹), and liquid saturation may remain far below 1.
Leaving these heterogeneous magnitudes directly in the solver vector can make
the numerical problem less well conditioned and can blur the practical meaning
of uniform solver tolerances.

The goal of the scaling is therefore purely numerical:
- bring most solver variables closer to order 1,
- improve conditioning of the time integration problem,
- make `reltol`/`atol` act more homogeneously across state components.

Important design choice
-----------------------
The physical model itself is not reformulated in scaled variables. The
electrochemical, transport and thermal equations remain written and evaluated in
physical units; only the solver-facing state vector is scaled/unscaled at the
numerical boundary.

The scaling is split by subsystem:
- cell: electrochemical core states (MEA + gas channels)
- manifold: manifold states (currently identity placeholders)
- auxiliary: BoP auxiliary states (currently identity placeholders)
- current_distribution: nonlinear GC current split solve (U_cell, i_fc, C_O2_Pt)
"""
module StateScalingModule

export CellStateScaling, ManifoldStateScaling, AuxiliaryStateScaling,
       CurrentDistributionScaling, DAEAlgebraicScaling, StateScaling

"""
Scaling references for cell state variables integrated by the solver.

These references are chosen so that the scaled cell states handled by the ODE
solver are usually of order 1, despite the heterogeneous physical magnitudes of
the original variables.
"""
Base.@kwdef struct CellStateScaling
    C_v::Float64 = 10.0
    C_H2::Float64 = 70.0
    C_O2::Float64 = 10.0
    C_N2::Float64 = 50.0
    T::Float64 = 350.0
    lambda::Float64 = 10.0
    eta_c::Float64 = 0.3
    phi::Float64 = 1.0
    C_CO2::Float64 = 10.0
    s::Float64 = 0.1
end

# Preserve the positional constructor exported before electrode potentials were
# added to the spatial state. The potential reference uses its default value.
CellStateScaling(C_v::Real, C_H2::Real, C_O2::Real, C_N2::Real, T::Real,
                 lambda::Real, eta_c::Real, s::Real) =
    CellStateScaling(Float64(C_v), Float64(C_H2), Float64(C_O2), Float64(C_N2),
                     Float64(T), Float64(lambda), Float64(eta_c), 1.0, 10.0, Float64(s))

"""
Scaling references for manifold states.
"""
Base.@kwdef struct ManifoldStateScaling
    P::Float64 = 1.0e5
    Phi::Float64 = 1.0
end

"""
Scaling references for auxiliary states.

Current values are identity placeholders. This keeps the auxiliary subsystem in
physical units until dedicated references are defined and validated.
"""
Base.@kwdef struct AuxiliaryStateScaling
    Wcp::Float64 = 1.0
    Wa_inj::Float64 = 1.0
    Wc_inj::Float64 = 1.0
    Abp_a::Float64 = 1.0
    Abp_c::Float64 = 1.0
end

"""
Scaling references for the GC current-distribution nonlinear solve.

These references are local to the algebraic system solved in
`calculate_1D_GC_current_density`, whose unknown vector is
`[U_cell, i_fc[1:nb_gc], C_O2_Pt[1:nb_gc]]`.
"""
Base.@kwdef struct CurrentDistributionScaling
    U::Float64 = 1.0
    i_fc::Float64 = 1.0e4
    C_O2_Pt::Float64 = 10.0
end

"""
Scaling references for algebraic variables appended in the IDA DAE state.

Unknown ordering is fixed as:
`[U_cell, i_fc[1:nb_gc], C_O2_Pt[1:nb_gc], J_a_in, J_c_in]`.
"""
Base.@kwdef struct DAEAlgebraicScaling
    U::Float64 = 1.0
    i_fc::Float64 = 1.0e4
    C_O2_Pt::Float64 = 10.0
    J_in::Float64 = 1.0e3
end

"""
Top-level state scaling configuration used internally by the solver interface.

It groups the fixed reference values applied to each subsystem when converting
between physical states and solver states.
"""
Base.@kwdef struct StateScaling
    cell::CellStateScaling = CellStateScaling()
    manifold::ManifoldStateScaling = ManifoldStateScaling()
    auxiliary::AuxiliaryStateScaling = AuxiliaryStateScaling()
    current_distribution::CurrentDistributionScaling = CurrentDistributionScaling()
    dae_algebraic::DAEAlgebraicScaling = DAEAlgebraicScaling()
end

CellStateScaling(cv,h2,o2,n2,T,lambda,eta,phi,s) =
    CellStateScaling(cv,h2,o2,n2,T,lambda,eta,phi,10.0,s)

end # module
