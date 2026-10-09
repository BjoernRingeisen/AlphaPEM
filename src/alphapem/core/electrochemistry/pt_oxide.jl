"""Finite sublayer inventory of spherical Pt particles (equivalent layers)."""
pt_bulk_capacity(p::ReactionParams) = max(p.Pt_particle_diameter / (4p.Pt_atomic_radius), 1.0)

"""One-layer Pt inventory in mol/m² of geometric electrode area."""
function pt_surface_inventory(side::Symbol, p::ReactionParams)
    side in (:acl, :ccl) || throw(ArgumentError("Pt inventory side must be acl or ccl"))
    loading = side == :acl ? p.Pt_loading_a : p.Pt_loading_c
    return 4loading * p.Pt_atomic_radius / (p.Pt_molar_mass * p.Pt_particle_diameter)
end

free_pt(c::PtCoverage) = 1 - c.OH - c.sO
function valid_pt_coverage(c::PtCoverage, p::ReactionParams; atol::Real=0.0)
    return all(isfinite, (c.OH, c.sO, c.bO)) &&
        c.OH >= -atol && c.sO >= -atol && free_pt(c) >= -atol &&
        -atol <= c.bO <= pt_bulk_capacity(p) + atol
end

# Extend the RHS to trial Newton states. Integrated states are never projected:
# conservation is carried by the reaction rates and a solver-domain guard.
function bounded_pt_coverage(c::PtCoverage, p::ReactionParams)
    all(isfinite, (c.OH, c.sO, c.bO)) || throw(ArgumentError("Pt coverage must be finite"))
    oh = clamp(c.OH, 0.0, 1.0)
    return PtCoverage(oh, clamp(c.sO, 0.0, 1 - oh), clamp(c.bO, 0.0, pt_bulk_capacity(p)))
end

"""Python's Pt → PtOH → surface PtO and gated bulk oxidation.

Oxidation currents are positive. Bulk oxidation does not consume surface PtO.
The surface balance closes exactly through the dependent free-Pt fraction.
"""
function pt_oxide(c::PtCoverage, phi::Real, T::Real, roughness::Real,
                  p::ReactionParams; side::Symbol=:ccl)
    validate_reaction_parameters(p)
    all(isfinite, (phi, T, roughness)) && T > 0 && roughness > 0 ||
        throw(ArgumentError("Invalid Pt-oxide reaction inputs"))
    c = bounded_pt_coverage(c, p)
    capacity = pt_bulk_capacity(p)
    remaining = max(capacity - c.bO, 0.0)
    capacity_gate = (remaining / (remaining + p.Pt_bO_capacity_smooth_layers))^p.Pt_bO_capacity_order
    surface_gate = 0.5 * (1 + tanh((c.sO - p.Pt_bO_surface_threshold) / p.Pt_bO_surface_width))
    interaction = safe_exp(p.Pt_bO_interaction_energy * c.bO / (R * T), p.exponent_limit)
    current(i0, Eact, E, alpha_ox, alpha_red, ox, red) = begin
        eta = F * (phi - E) / (R * T)
        roughness * i0 * safe_exp(Eact / R * (1 / p.Tref_PtOx - 1 / T), 80.0) *
            (ox * safe_exp(alpha_ox * eta, p.exponent_limit) -
             red * safe_exp(-alpha_red * eta, p.exponent_limit))
    end
    j1 = current(p.i0_PtOH, p.Eact_PtOH, p.E_PtOH,
                 p.alpha_z_PtOH_ox, p.alpha_z_PtOH_red, free_pt(c), c.OH)
    j2 = current(p.i0_Pt_sO, p.Eact_Pt_sO, p.E_Pt_sO,
                 p.alpha_z_Pt_sO_ox, p.alpha_z_Pt_sO_red, c.OH, c.sO)
    j3 = current(p.i0_Pt_bO, p.Eact_Pt_bO, p.E_Pt_bO,
                 p.alpha_z_Pt_bO_ox, p.alpha_z_Pt_bO_red,
                 c.sO * surface_gate * capacity_gate * interaction, c.bO)
    inventory = pt_surface_inventory(side, p)
    charge = F * inventory
    return PtOxideResult((j1, j2, j3), PtCoverage((j1 - j2) / charge,
                         j2 / charge, j3 / (2charge)), inventory, capacity)
end

with_potential(s::ElectrodeState, phi::Real) =
    ElectrodeState(s.T, s.C_H2, s.C_O2, phi, s.thickness, s.coverage,s.C_CO2,s.enable_cor)

"""Stationary oxide coverages at fixed potential with finite bulk capacity."""
function stationary_pt_coverage(phi::Real, T::Real, p::ReactionParams;
        C_CO2::Real=0.0,enable_cor::Bool=false,side::Symbol=:ccl)
    validate_reaction_parameters(p)
    isfinite(phi) && isfinite(T) && T > 0 || throw(ArgumentError("Invalid stationary oxide inputs"))
    exponent(alpha, E) = clamp(alpha * F * (phi - E) / (R * T), -p.exponent_limit, p.exponent_limit)
    l1 = exponent(p.alpha_z_PtOH_ox, p.E_PtOH) - exponent(-p.alpha_z_PtOH_red, p.E_PtOH)
    l2 = exponent(p.alpha_z_Pt_sO_ox, p.E_Pt_sO) - exponent(-p.alpha_z_Pt_sO_red, p.E_Pt_sO)
    if enable_cor
        p.i0_PtOH == 0 && return PtCoverage()
        unit = ElectrodeState(T,0.0,0.0,phi,1.0,PtCoverage(1,0,0),C_CO2,true)
        cor = carbon_oxidation(unit,p;side)
        roughness = side == :acl ? p.roughness_a : p.roughness_c
        k1 = roughness*p.i0_PtOH*safe_exp(p.Eact_PtOH/R*(1/p.Tref_PtOx-1/T),80.0)
        # j_PtOH = j_COR_cat/3 when surface PtO is stationary.
        lred = exponent(-p.alpha_z_PtOH_red,p.E_PtOH)
        lcat = cor.cat.current == 0 ? -Inf : log(cor.cat.current/(3k1))
        m = max(lred,lcat)
        l1 = exponent(p.alpha_z_PtOH_ox,p.E_PtOH) - (m+log(exp(lred-m)+exp(lcat-m)))
    end
    weights = (0.0, l1, l1 + l2)
    scale = maximum(weights)
    w = exp.(weights .- scale)
    oh, so = w[2] / sum(w), w[3] / sum(w)
    lo, hi = 0.0, pt_bulk_capacity(p)
    for _ in 1:100
        mid = (lo + hi) / 2
        j = pt_oxide(PtCoverage(oh, so, mid), phi, T, 1.0, p).currents[3]
        j == 0 && return PtCoverage(oh, so, mid)
        if j > 0
            lo = mid
        else
            hi = mid
        end
    end
    return PtCoverage(oh, so, (lo + hi) / 2)
end

"""Initialize potential with prescribed or jointly stationary Pt coverages.

Prescribed clean Pt follows Python's default explicit-startup preparation.
Stationary mode solves Faraday balance and all coverage balances together.
"""
function initialize_pt_electrode(s::ElectrodeState, p::ReactionParams;
        side::Symbol=:ccl, mode::Symbol=:prescribed_coverages,
        coverage::PtCoverage=PtCoverage(), external_current::Real=0.0,
        C_O2_Pt::Real=s.C_O2, transport_ratio::Union{Nothing,Real}=nothing)
    mode in (:prescribed_coverages, :stationary_local) || throw(ArgumentError("Invalid Pt initialization mode"))
    valid_pt_coverage(coverage, p) || throw(ArgumentError("Invalid prescribed Pt coverage"))
    roughness = side == :acl ? p.roughness_a : p.roughness_c
    initial = ElectrodeState(s.T, s.C_H2, s.C_O2, s.phi, s.thickness, coverage,s.C_CO2,s.enable_cor)
    phi = stationary_potential(initial, roughness, p; external_current, C_O2_Pt, side, transport_ratio,
                                stationary_coverages=mode == :stationary_local)
    cov = mode == :stationary_local ? stationary_pt_coverage(phi, s.T, p;
        C_CO2=s.C_CO2,enable_cor=s.enable_cor,side) : coverage
    return ElectrodeState(s.T, s.C_H2, s.C_O2, phi, s.thickness, cov,s.C_CO2,s.enable_cor)
end
