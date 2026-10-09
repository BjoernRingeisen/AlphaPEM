"""Local HOR/ORR/Pt-oxide kinetics and conservative catalyst-layer balances.

HOR currents are oxidation-positive; ORR currents are reduction-positive.
The net Faraday current is HOR minus ORR plus enabled Pt-oxide currents. All returned
currents use geometric area. No plotting, solver or Python dependency is needed.
"""
module Electrochemistry

using ...Config: ReactionParams, validate_reaction_parameters
using ...Utils: F, R, E0, Pref_eq, delta_s_HOR, delta_s_ORR
using ..Types: ElectrodeState, ReactionResult, ElectrodeReactions, MultireactionResult, PtCoverage, PtOxideResult
using ..Types: CarbonOxidationResult

export ElectrodeState, ReactionResult, ElectrodeReactions, MultireactionResult,
       reversible_hor_potential, reversible_orr_potential, hor_current, orr_current,
       electrode_reactions, multireaction, reaction_sources, reaction_heat,
       potential_derivatives, stationary_potential, limit_inventory_current,
       pt_oxygen_residual, resolve_pt_oxygen, electrode_dae_residual!,
       PtCoverage, PtOxideResult, pt_oxide, pt_bulk_capacity, pt_surface_inventory,
       valid_pt_coverage, free_pt, stationary_pt_coverage, initialize_pt_electrode
export CarbonOxidationResult, carbon_oxidation, reversible_cor_potential

include("pt_oxide.jl")
include("carbon_oxidation.jl")

@inline safe_exp(x, limit) = exp(clamp(x, -limit, limit))
@inline positive(c) = max(c, 0.0)
@inline availability(c, p) = (positive(c) / (positive(c) + p.concentration_scale))^p.availability_order
@inline gas_pressure(c, T, p) = max(R * T * positive(c), R * T * p.concentration_floor, 1e-12)

function check_state(s::ElectrodeState)
    all(isfinite, (s.T, s.C_H2, s.C_O2, s.C_CO2, s.phi, s.thickness)) ||
        throw(ArgumentError("Electrode state must be finite"))
    s.T > 0 && s.thickness > 0 || throw(ArgumentError("Temperature and CL thickness must be positive"))
    return s
end

function reversible_hor_potential(c::Real, T::Real, p::ReactionParams)
    validate_reaction_parameters(p)
    T > 0 && isfinite(T) && isfinite(c) || throw(ArgumentError("Invalid Nernst inputs"))
    return -R * T / (2F) * log(gas_pressure(c, T, p) / Pref_eq)
end

function reversible_orr_potential(c::Real, T::Real, p::ReactionParams)
    validate_reaction_parameters(p)
    T > 0 && isfinite(T) && isfinite(c) || throw(ArgumentError("Invalid Nernst inputs"))
    return E0 - 8.5e-4 * (T - 298.15) + R * T / (4F) * log(gas_pressure(c, T, p) / Pref_eq)
end

"""Smooth cap on gas-consuming current; reverse (gas-producing) current is unchanged."""
function limit_inventory_current(j::Real, c::Real, H::Real, n::Real, p::ReactionParams)
    validate_reaction_parameters(p)
    all(isfinite, (j, c, H, n)) && H > 0 && n > 0 || throw(ArgumentError("Invalid inventory limiter inputs"))
    j <= 0 && return (current=Float64(j), factor=1.0)
    cap = positive(c) * n * F * H / p.inventory_time
    cap <= 1e-30 && return (current=0.0, factor=0.0)
    # Equivalent to j/(1+(j/cap)^order)^(1/order), without overflowing j/cap.
    low, high = min(j, cap), max(j, cap)
    limited = low / (1 + (low / high)^p.inventory_order)^(1 / p.inventory_order)
    return (current=limited, factor=limited / j)
end

function kinetic_current(c, eta, T, roughness, i0, Eact, Tref, cref, alpha_forward, alpha_reverse, p)
    validate_reaction_parameters(p)
    all(isfinite, (c, eta, T, roughness)) && T > 0 && roughness > 0 ||
        throw(ArgumentError("Invalid kinetic inputs"))
    cpos = positive(c)
    cpos == 0 && return 0.0
    temperature = safe_exp(Eact / R * (1 / Tref - 1 / T), 80.0)
    x = F * eta / (R * T)
    bracket = safe_exp(alpha_forward * x, p.exponent_limit) - safe_exp(-alpha_reverse * x, p.exponent_limit)
    return availability(c, p) * (cpos / cref)^p.concentration_power * roughness * i0 * temperature * bracket
end

"""HOR current before the inventory cap, with oxidation overpotential phi-E_HOR."""
hor_current(c, eta, T, roughness, p::ReactionParams) = kinetic_current(
    c, eta, T, roughness, p.i0_hor, p.Eact_hor, p.Tref_hor, p.C_H2_ref,
    p.alpha_z_hor_ox, p.alpha_z_hor_red, p)

"""ORR current before the inventory cap, with reduction overpotential E_ORR-phi."""
orr_current(c, eta, T, roughness, p::ReactionParams) = kinetic_current(
    c, eta, T, roughness, p.i0_orr, p.Eact_orr, p.Tref_orr, p.C_O2_ref,
    p.alpha_z_orr_red, p.alpha_z_orr_ox, p)

function electrode_reactions(s::ElectrodeState, roughness::Real, p::ReactionParams;
                             C_O2_Pt::Real=s.C_O2, side::Symbol=:ccl)
    validate_reaction_parameters(p)
    check_state(s)
    E_HOR = reversible_hor_potential(s.C_H2, s.T, p)
    E_ORR = reversible_orr_potential(C_O2_Pt, s.T, p)
    eta_HOR, eta_ORR = s.phi - E_HOR, E_ORR - s.phi
    raw_hor = hor_current(s.C_H2, eta_HOR, s.T, roughness, p)
    raw_orr = orr_current(C_O2_Pt, eta_ORR, s.T, roughness, p)
    oxide = s.coverage === nothing ? nothing : pt_oxide(s.coverage, s.phi, s.T, roughness, p; side)
    s.enable_cor && s.coverage === nothing && throw(ArgumentError("COR requires Pt coverage states"))
    cor = s.enable_cor ? carbon_oxidation(s,p;side) : nothing
    if cor !== nothing
        rates = oxide.rates
        oxide = PtOxideResult(oxide.currents,
            PtCoverage(rates.OH-cor.ptOH_consumption_rate,rates.sO,rates.bO),
            oxide.surface_inventory,oxide.max_layers)
    end
    if s.coverage !== nothing
        coverage = bounded_pt_coverage(s.coverage, p)
        # Python treats Pt and PtOH as equally ORR-active, and PtO as a
        # separate, less active branch. HOR retains its source-model law.
        raw_orr = kinetic_current(C_O2_Pt, eta_ORR, s.T, roughness,
            (1 - coverage.sO) * p.i0_orr + coverage.sO * p.i0_orr_pto,
            p.Eact_orr, p.Tref_orr, p.C_O2_ref, p.alpha_z_orr_red, p.alpha_z_orr_ox, p)
    end
    hor = limit_inventory_current(raw_hor, s.C_H2, s.thickness, 2, p)
    # Oxygen is removed from the bulk CL; Pt concentration controls kinetics only.
    orr = limit_inventory_current(raw_orr, s.C_O2, s.thickness, 4, p)
    return ElectrodeReactions(
        ReactionResult(hor.current, raw_hor, E_HOR, eta_HOR, availability(s.C_H2, p), hor.factor, s.C_H2),
        ReactionResult(orr.current, raw_orr, E_ORR, eta_ORR, availability(C_O2_Pt, p), orr.factor, C_O2_Pt),
        hor.current - orr.current + (oxide === nothing ? 0.0 : sum(oxide.currents)) +
            (cor === nothing ? 0.0 : cor.noncat.current+cor.cat.current), oxide, cor)
end

"""Evaluate both electrodes using a caller-supplied, transport-consistent Pt concentration."""
function multireaction(a::ElectrodeState, c::ElectrodeState, p::ReactionParams; C_O2_Pt::Real=c.C_O2)
    validate_reaction_parameters(p)
    return MultireactionResult(electrode_reactions(a, p.roughness_a, p; side=:acl),
                              electrode_reactions(c, p.roughness_c, p; C_O2_Pt, side=:ccl))
end

"""Double-layer derivatives (V/s); external local current is +i at ACL and -i at CCL."""
function potential_derivatives(r::MultireactionResult, i::Real, p::ReactionParams)
    validate_reaction_parameters(p)
    isfinite(i) || throw(ArgumentError("Current must be finite"))
    return (phi_a=(i - r.anode.faraday) / p.Cdl_a,
            phi_c=(-i - r.cathode.faraday) / p.Cdl_c)
end

"""Reaction and crossover sources in mol/m³/s.

J_H2 and J_O2 are signed membrane fluxes in mol/m²/s, positive ACL to CCL.
They transfer gas; reactions consume it locally. Water enters the ionomer
balance. Do not add the legacy instantaneous crossover reaction sources.
"""
function reaction_sources(r::MultireactionResult, H_a::Real, H_c::Real; J_H2::Real=0.0, J_O2::Real=0.0)
    all(isfinite, (H_a, H_c, J_H2, J_O2)) && H_a > 0 && H_c > 0 ||
        throw(ArgumentError("Invalid source geometry or flux"))
    sources(e, H, sign) = (
        H2=-e.hor.current / (2F * H) + sign * J_H2 / H,
        O2=-e.orr.current / (4F * H) + sign * J_O2 / H,
        H2O=e.orr.current / (2F * H) - (e.oxide === nothing ? 0.0 :
            (e.oxide.currents[1] + e.oxide.currents[3] / 2) / (F * H)) -
            (e.cor === nothing ? 0.0 : e.cor.noncat.current/(2F*H)+e.cor.cat.current/(3F*H)),
        CO2=e.cor === nothing ? 0.0 : e.cor.noncat.current/(4F*H)+e.cor.cat.current/(3F*H),
        carbon=e.cor === nothing ? 0.0 : -(e.cor.noncat.current/(4F*H)+e.cor.cat.current/(3F*H)))
    return (anode=sources(r.anode, H_a, -1), cathode=sources(r.cathode, H_c, 1))
end

"""Volumetric entropy, gas-activity and activation heat for one CL (W/m³).

Uses the source model's liquid-water activity of one and proton reference.
Joule and sorption heat remain separate transport-model contributions. Direct
Pt-oxide heat is excluded, matching Python's unspecified oxide thermochemistry;
the coverage-dependent gas reaction currents still enter all heat terms.
The gas activity defaults to the concentration recorded by the reaction
evaluation, not bulk cathode oxygen. Passing a different Pt concentration
is an error: re-evaluate the reactions for the new state instead.
"""
function reaction_heat(s::ElectrodeState, e::ElectrodeReactions, p::ReactionParams; C_O2_Pt::Real=e.orr.concentration)
    validate_reaction_parameters(p)
    check_state(s)
    isfinite(C_O2_Pt) && C_O2_Pt == e.orr.concentration ||
        throw(ArgumentError("Heat must use the Pt concentration of the evaluated reaction"))
    s.C_H2 == e.hor.concentration || throw(ArgumentError("Hydrogen state differs from evaluated reaction"))
    jh, jo = e.hor.current, e.orr.current
    entropy = s.T / (F * s.thickness) * (-jh * delta_s_HOR / 2 - jo * delta_s_ORR / 4)
    logH2 = -log(gas_pressure(s.C_H2, s.T, p) / Pref_eq)
    logO2 = -log(gas_pressure(C_O2_Pt, s.T, p) / Pref_eq)
    nernst = R * s.T / (F * s.thickness) * (jh * logH2 / 2 + jo * logO2 / 4)
    activation = (jh * e.hor.overpotential + jo * e.orr.overpotential) / s.thickness
    if e.cor !== nothing
        s.C_CO2 == e.cor.noncat.concentration || throw(ArgumentError("CO₂ state differs from evaluated COR"))
        for (r,n,ds) in ((e.cor.noncat,4,p.delta_s_cor_noncat),(e.cor.cat,3,p.delta_s_cor_cat))
            entropy += -r.current*s.T*ds/(n*F*s.thickness)
            nernst += r.current*R*s.T*log(gas_pressure(s.C_CO2,s.T,p)/Pref_eq)/(n*F*s.thickness)
            activation += r.current*r.overpotential/s.thickness
        end
    end
    return (entropy=entropy, nernst=nernst, activation=activation, total=entropy + nernst + activation)
end

"""Cathode Pt-oxygen algebraic residual in mol/m³.

`transport_ratio` is R_T_O2_Pt/a_c in seconds, supplied by the host transport
model. The ORR Faraday current sets oxygen consumption, independently of the
external current. The floor reproduces the Python concentration constraint.
"""
function pt_oxygen_residual(C_O2_Pt::Real, s::ElectrodeState, p::ReactionParams,
                            transport_ratio::Real; concentration_floor::Real=1e-6)
    all(isfinite, (transport_ratio, concentration_floor)) && transport_ratio >= 0 && concentration_floor > 0 ||
        throw(ArgumentError("Invalid Pt oxygen transport ratio or floor"))
    r = electrode_reactions(s, p.roughness_c, p; C_O2_Pt)
    bulk = max(s.C_O2, concentration_floor)
    target = max(bulk - max(r.orr.current, 0.0) * transport_ratio / (4F * s.thickness), concentration_floor)
    return C_O2_Pt - target
end

"""Solve the local Pt-oxygen constraint by a bracketed solve.

Uses the Python transport law but checks convergence instead of accepting a
fixed number of relaxed iterations. This is useful for local initialization;
the spatial DAE can use `pt_oxygen_residual` directly.
"""
function resolve_pt_oxygen(s::ElectrodeState, p::ReactionParams, transport_ratio::Real;
                          concentration_floor::Real=1e-6, atol::Real=1e-10, maxiters::Int=100)
    isfinite(atol) && atol > 0 && maxiters > 0 || throw(ArgumentError("Invalid Pt solve controls"))
    residual(c) = pt_oxygen_residual(c, s, p, transport_ratio; concentration_floor)
    lo, hi = Float64(concentration_floor), max(s.C_O2, concentration_floor)
    flo, fhi = residual(lo), residual(hi)
    abs(flo) <= atol && return lo
    abs(fhi) <= atol && return hi
    flo < 0 < fhi || throw(ArgumentError("No Pt oxygen concentration bracket"))
    for _ in 1:maxiters
        mid = (lo + hi) / 2
        fmid = residual(mid)
        abs(fmid) <= atol && return mid
        if fmid < 0
            lo = mid
        else
            hi = mid
        end
    end
    error("Pt oxygen concentration solve did not converge")
end

"""Write the local physical-unit DAE residual for [phi_a, phi_c, C_O2_Pt].

The first two states are differential; Pt oxygen is algebraic. The supplied
electrode states provide T, gas concentrations and CL thickness; their phi
fields are replaced by `y`. `dy[3]` is deliberately unused. This adapter does
not add spatial inventories or enable a new `run_simulation` mode.
"""
function electrode_dae_residual!(res::AbstractVector, dy::AbstractVector, y::AbstractVector,
                                  a::ElectrodeState, c::ElectrodeState, i::Real,
                                  p::ReactionParams, transport_ratio::Real)
    length(res) == length(dy) == length(y) == 3 || throw(DimensionMismatch("Expected three local DAE states"))
    anode = with_potential(a, y[1])
    cathode = with_potential(c, y[2])
    r = multireaction(anode, cathode, p; C_O2_Pt=y[3])
    d = potential_derivatives(r, i, p)
    res[1] = dy[1] - d.phi_a
    res[2] = dy[2] - d.phi_c
    res[3] = pt_oxygen_residual(y[3], cathode, p, transport_ratio)
    return nothing
end

"""Solve local Faraday balance at fixed inventories by bounded bisection.

The external current is oxidation-positive (use -i at the cathode). A target
outside the achievable interval raises an error rather than inventing an
initial potential. Optional stationary coverages and Pt transport are solved
inside the potential residual for coupled initialization.
"""
function stationary_potential(s::ElectrodeState, roughness::Real, p::ReactionParams;
                              external_current::Real=0.0, C_O2_Pt::Real=s.C_O2,
                              bounds=(-0.2, 1.5), atol::Real=1e-8, maxiters::Int=100,
                              side::Symbol=:ccl, stationary_coverages::Bool=false,
                              transport_ratio::Union{Nothing,Real}=nothing)
    validate_reaction_parameters(p)
    check_state(s)
    lo, hi = Float64.(bounds)
    all(isfinite, (lo, hi, external_current, atol)) && lo < hi && atol > 0 && maxiters > 0 ||
        throw(ArgumentError("Invalid stationary-potential solve controls"))
    residual(phi) = begin
        coverage = stationary_coverages ? stationary_pt_coverage(phi, s.T, p;
            C_CO2=s.C_CO2,enable_cor=s.enable_cor,side) : s.coverage
        state = ElectrodeState(s.T, s.C_H2, s.C_O2, phi, s.thickness, coverage,s.C_CO2,s.enable_cor)
        pt = transport_ratio === nothing ? C_O2_Pt : resolve_pt_oxygen(state, p, transport_ratio)
        electrode_reactions(state, roughness, p; C_O2_Pt=pt, side).faraday - external_current
    end
    flo, fhi = residual(lo), residual(hi)
    abs(flo) <= atol && abs(fhi) <= atol && throw(ArgumentError("Potential is undetermined at these inventories"))
    abs(flo) <= atol && return lo
    abs(fhi) <= atol && return hi
    flo < 0 < fhi || throw(ArgumentError("No stationary potential bracket for current $external_current"))
    for _ in 1:maxiters
        mid = (lo + hi) / 2
        fmid = residual(mid)
        abs(fmid) <= atol && return mid
        if fmid < 0
            lo = mid
        else
            hi = mid
        end
    end
    error("Stationary-potential solve did not converge")
end

end
