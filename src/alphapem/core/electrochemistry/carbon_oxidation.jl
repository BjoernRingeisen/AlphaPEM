"""Reduction potential of COR vs SHE, including the CO₂ gas activity."""
function reversible_cor_potential(c::Real, T::Real, p::ReactionParams; catalyzed::Bool=false)
    validate_reaction_parameters(p)
    T > 0 && isfinite(T) && isfinite(c) || throw(ArgumentError("Invalid COR Nernst inputs"))
    n = catalyzed ? 3.0 : 4.0
    dh = catalyzed ? p.delta_h_cor_cat : p.delta_h_cor_noncat
    ds = catalyzed ? p.delta_s_cor_cat : p.delta_s_cor_noncat
    return -(dh - T*ds)/(n*F) - p.SHE_absolute_potential +
           R*T/(n*F)*log(gas_pressure(c,T,p)/Pref_eq)
end

"""Python's irreversible 4-electron carbon and 3-electron PtOH-catalyzed COR.

The catalyzed reaction consumes one PtOH site per carbon atom, returning free
Pt. The source model assumes available solid carbon and unit water activity.
"""
function carbon_oxidation(s::ElectrodeState, p::ReactionParams; side::Symbol=:ccl)
    validate_reaction_parameters(p)
    check_state(s)
    side in (:acl,:ccl) || throw(ArgumentError("Invalid COR electrode side"))
    oh = s.coverage === nothing ? 0.0 : bounded_pt_coverage(s.coverage,p).OH
    rc = side == :acl ? p.roughness_carbon_a : p.roughness_carbon_c
    rp = side == :acl ? p.roughness_a : p.roughness_c
    reaction(cat) = begin
        eq = reversible_cor_potential(s.C_CO2,s.T,p; catalyzed=cat)
        eta = s.phi - eq
        i0 = cat ? p.i0_cor_cat : p.i0_cor_noncat
        ea = cat ? p.Eact_cor_cat : p.Eact_cor_noncat
        alpha = cat ? p.alpha_z_cor_cat : p.alpha_z_cor_noncat
        activity = (cat ? oh : 1.0)*p.cor_activity_factor
        raw = (cat ? rp : rc)*(cat ? oh : 1.0)*i0*
            safe_exp(ea/R*(1/p.Tref_cor-1/s.T),80.0)*
            safe_exp(alpha*F*eta/(R*s.T),p.exponent_limit)
        ReactionResult(p.cor_activity_factor*raw,raw,eq,eta,activity,1.0,s.C_CO2)
    end
    noncat, cat = reaction(false), reaction(true)
    carbon_rate = noncat.current/(4F) + cat.current/(3F)
    return CarbonOxidationResult(noncat,cat,p.carbon_molar_mass*carbon_rate,
        cat.current/(3F*pt_surface_inventory(side,p)))
end
