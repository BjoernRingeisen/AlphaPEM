"""HOR/ORR parameters ported from Python's multireaction model at d7a65bb.

Exchange currents are A/m² of Pt; roughness converts them to geometric area.
Capacitances are F/m² of geometric area. Concentrations are mol/m³.
These source-branch defaults are a reproduction set, not a fit to every stack.
"""
Base.@kwdef struct ReactionParams
    i0_hor::Float64 = 2000.0
    i0_orr::Float64 = 7.3e-4
    Eact_hor::Float64 = 10000.0
    Eact_orr::Float64 = 55000.0
    Tref_hor::Float64 = 353.15
    Tref_orr::Float64 = 353.15
    alpha_z_hor_ox::Float64 = 0.3
    alpha_z_hor_red::Float64 = 0.3
    alpha_z_orr_ox::Float64 = 0.65
    alpha_z_orr_red::Float64 = 1.0
    C_H2_ref::Float64 = 39.195263103187315
    C_O2_ref::Float64 = 6.194215337579467
    concentration_power::Float64 = 2.0
    concentration_floor::Float64 = 1e-30
    concentration_scale::Float64 = 1e-6
    availability_order::Float64 = 4.0
    inventory_time::Float64 = 1e-2
    inventory_order::Float64 = 2.0
    exponent_limit::Float64 = 80.0
    roughness_a::Float64 = 66.0
    roughness_c::Float64 = 200.0
    Cdl_a::Float64 = 130.0
    Cdl_c::Float64 = 400.0
    i0_orr_pto::Float64 = 3.5e-4
    E_PtOH::Float64 = 0.730
    E_Pt_sO::Float64 = 0.850
    E_Pt_bO::Float64 = 0.850
    i0_PtOH::Float64 = 6.0
    i0_Pt_sO::Float64 = 8.0
    i0_Pt_bO::Float64 = 13.0
    Eact_PtOH::Float64 = 2670.0
    Eact_Pt_sO::Float64 = 94830.0
    Eact_Pt_bO::Float64 = 94830.0
    Tref_PtOx::Float64 = 353.15
    alpha_z_PtOH_ox::Float64 = 0.45
    alpha_z_PtOH_red::Float64 = 0.30
    alpha_z_Pt_sO_ox::Float64 = 0.20
    alpha_z_Pt_sO_red::Float64 = 0.20
    alpha_z_Pt_bO_ox::Float64 = 0.60
    alpha_z_Pt_bO_red::Float64 = 0.40
    Pt_bO_interaction_energy::Float64 = 1900.0
    Pt_bO_surface_threshold::Float64 = 0.0
    Pt_bO_surface_width::Float64 = 5e-3
    Pt_bO_capacity_smooth_layers::Float64 = 0.05
    Pt_bO_capacity_order::Float64 = 2.0
    Pt_loading_a::Float64 = 0.002
    Pt_loading_c::Float64 = 0.006
    Pt_particle_diameter::Float64 = 4e-9
    Pt_atomic_radius::Float64 = 0.139e-9
    Pt_molar_mass::Float64 = 0.195084
    roughness_carbon_a::Float64 = 264.0
    roughness_carbon_c::Float64 = 800.0
    Tref_cor::Float64 = 353.15
    SHE_absolute_potential::Float64 = 14.66
    i0_cor_noncat::Float64 = 1e-19
    i0_cor_cat::Float64 = 1e-7
    Eact_cor_noncat::Float64 = 76700.0
    Eact_cor_cat::Float64 = 76700.0
    alpha_z_cor_noncat::Float64 = 0.90
    alpha_z_cor_cat::Float64 = 0.52
    delta_h_cor_noncat::Float64 = -5837890.0
    delta_h_cor_cat::Float64 = -4304790.0
    delta_s_cor_noncat::Float64 = -329.26
    delta_s_cor_cat::Float64 = 0.0
    carbon_molar_mass::Float64 = 0.012011
    cor_activity_factor::Float64 = 1.0
end

function validate_reaction_parameters(p::ReactionParams)
    nonnegative = (:i0_hor, :i0_orr, :Eact_hor, :Eact_orr, :concentration_power,
                   :i0_orr_pto, :i0_PtOH, :i0_Pt_sO, :i0_Pt_bO,
                   :Eact_PtOH, :Eact_Pt_sO, :Eact_Pt_bO, :Pt_bO_surface_threshold,
                   :i0_cor_noncat, :i0_cor_cat, :Eact_cor_noncat, :Eact_cor_cat,
                   :cor_activity_factor)
    signed = (:delta_h_cor_noncat, :delta_h_cor_cat, :delta_s_cor_noncat, :delta_s_cor_cat)
    for name in fieldnames(ReactionParams)
        value = getfield(p, name)
        valid = isfinite(value) && (name in signed || (name in nonnegative ? value >= 0 : value > 0))
        valid || throw(ArgumentError("Invalid reaction parameter $name: $value"))
    end
    p.exponent_limit <= 80 || throw(ArgumentError("exponent_limit must not exceed 80"))
    p.Pt_bO_surface_threshold <= 1 || throw(ArgumentError("PtO threshold must be at most one"))
    p.Pt_bO_capacity_order >= 1 || throw(ArgumentError("Bulk capacity order must be at least one"))
    p.cor_activity_factor <= 1 || throw(ArgumentError("COR activity factor must be in [0, 1]"))
    return p
end
