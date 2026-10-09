# -*- coding: utf-8 -*-

"""This file represents all the differential equations used for the fuel cell model.
"""


# ____________________________________________________Main functions____________________________________________________

"""Convert a molar balance into a concentration rate for a changing gas-filled pore volume.

For inventory `epsilon * (1 - s) * C`, this includes the concentration change
caused by liquid saturation changing at rate `ds_dt`. Structural porosity is
allowed to change through `d_epsilon_dt`, which is nonzero in catalyst layers
when ionomer hydration or temperature changes their porosity.
"""
@inline function gas_concentration_rate(C, rhs, epsilon, s, ds_dt, d_epsilon_dt=0.0)
    gas_fraction = 1 - s
    return rhs / (epsilon * gas_fraction) + C * ds_dt / gas_fraction -
           C * d_epsilon_dt / epsilon
end

"""Directional time derivative of catalyst-layer porosity from lambda and temperature rates."""
@inline function cl_porosity_rate(element, lambda_cl, T_cl, Hcl, pp, dlambda_dt, dT_dt)
    (iszero(dlambda_dt) && iszero(dT_dt)) && return 0.0
    delta_t = 1e-6 / max(abs(dlambda_dt), abs(dT_dt), 1.0)
    epsilon_plus = epsilon_cl(element, lambda_cl + delta_t * dlambda_dt,
                              T_cl + delta_t * dT_dt, Hcl, pp)
    epsilon_minus = epsilon_cl(element, lambda_cl - delta_t * dlambda_dt,
                               T_cl - delta_t * dT_dt, Hcl, pp)
    return (epsilon_plus - epsilon_minus) / (2delta_t)
end

"""Temperature derivative of the liquid-density correlation, including its clamp."""
@inline function liquid_density_temperature_slope(T)
    T==Utils._liquid_water_temperature_value(T) || return 0.0
    c=T-273.15
    numerator_slope=16.945176-2*7.9870401e-3*c-3*46.170461e-6*c^2+
        4*105.56302e-9*c^3-5*280.54253e-12*c^4
    return (numerator_slope-rho_H2O_l(T)*16.879850e-3)/(1+16.879850e-3*c)
end

"""Conserve liquid mass epsilon*s*rho(T) as pore volume and temperature change."""
@inline function liquid_saturation_rate(s, mass_rhs, epsilon, T, dT_dt=0.0, d_epsilon_dt=0.0)
    rho=rho_H2O_l(T)
    return mass_rhs/(rho*epsilon)-s*d_epsilon_dt/epsilon-
        s*liquid_density_temperature_slope(T)*dT_dt/rho
end

"""Calculate dissolved-water (lambda) dynamics contribution.

Parameters
----------
sv : CellState1D{NB_GDL, NB_MPL}
    Typed 1D internal state for one gas-channel column.
pp : PhysicalParams
    Fuel-cell physical parameters container (geometry, thicknesses and porous properties).
S_abs : MEASorptionSources
    Water absorption/desorption rates at the CL ionomer (mol·m⁻³·s⁻¹).
J_lambda : MEADissolvedWaterFlux
    Dissolved-water inter-layer fluxes (mol·m⁻²·s⁻¹).
Sp : MEAWaterProductionSources
    Water production rates at the CLs (mol·m⁻³·s⁻

Returns
-------
MEADissolvedWaterDerivative
    Container with the lambda derivatives for the ACL, membrane and CCL.
"""
function calculate_dyn_dissoved_water_evolution_inside_MEA(
        sv::CellState1D{NB_GDL, NB_MPL},
        pp::PhysicalParams,
        S_abs::MEASorptionSources,
        J_lambda::MEADissolvedWaterFlux,
        Sp::MEAWaterProductionSources
)::MEADissolvedWaterDerivative where {NB_GDL, NB_MPL}

    # Extraction of the variables
    T_acl, T_ccl = sv.acl.T, sv.ccl.T
    lambda_acl, lambda_mem, lambda_ccl = sv.acl.lambda, sv.mem.lambda, sv.ccl.lambda
    M_eq, rho_mem = pp.M_eq, pp.rho_mem

    # Differential equations
    rhs_lambda_acl = -J_lambda.acl_mem / pp.Hacl + S_abs.v_acl + S_abs.l_acl + Sp.acl
    rhs_lambda_mem = (J_lambda.acl_mem - J_lambda.mem_ccl) / pp.Hmem
    rhs_lambda_ccl = J_lambda.mem_ccl / pp.Hccl + S_abs.v_ccl + S_abs.l_ccl + Sp.ccl

    C_fix_acl = cl_dry_ionomer_storage_capacity(:acl, pp.Hacl, pp)
    C_fix_ccl = cl_dry_ionomer_storage_capacity(:ccl, pp.Hccl, pp)

    d_lambda_acl_dt = rhs_lambda_acl / C_fix_acl
    d_lambda_mem_dt = M_eq / rho_mem * rhs_lambda_mem
    d_lambda_ccl_dt = rhs_lambda_ccl / C_fix_ccl

    return MEADissolvedWaterDerivative(d_lambda_acl_dt, d_lambda_mem_dt, d_lambda_ccl_dt)
end


"""Calculate the dynamic evolution of liquid water in the porous layers.

Parameters
----------
sv : CellState1D{NB_GDL, NB_MPL}
    Typed 1D internal state for one gas-channel column.
pp : PhysicalParams
    Fuel-cell physical parameters container (geometry, thicknesses and porous properties).
Jl : MEALiquidFluxes{NB_GDL, NB_MPL}
    Liquid-water inter-layer fluxes (kg·m⁻²·s⁻¹).
S_abs : MEASorptionSources
    Water absorption/desorption rates at the CL ionomer (mol·m⁻³·s⁻¹).
Sl : MEALiquidSources{NB_GDL, NB_MPL}
    Liquid-water phase-change source terms (mol·m⁻³·s⁻¹).

Returns
-------
CellDerivative1D{NB_GDL, NB_MPL}
    Updated derivative container with s (liquid saturation) derivatives filled in.
"""
function calculate_dyn_liquid_water_evolution_inside_MEA(
        sv::CellState1D{NB_GDL, NB_MPL},
        pp::PhysicalParams,
        Jl::MEALiquidFluxes{NB_GDL, NB_MPL},
        S_abs::MEASorptionSources,
        Sl::MEALiquidSources{NB_GDL, NB_MPL};
        temperature_derivative::Union{Nothing,MEATemperatureDerivative{NB_GDL,NB_MPL}}=nothing,
        cl_porosity_derivative=(acl=0.0,ccl=0.0)
)::MEALiquidWaterDerivative{NB_GDL, NB_MPL} where {NB_GDL, NB_MPL}

    # Extraction of the variables
    T_agdl, T_ampl = getproperty.(sv.agdl, :T), getproperty.(sv.ampl, :T)
    T_cmpl, T_cgdl = getproperty.(sv.cmpl, :T), getproperty.(sv.cgdl, :T)
    T_acl, T_ccl = sv.acl.T, sv.ccl.T
    lambda_acl, lambda_ccl = sv.acl.lambda, sv.ccl.lambda

    # Surface-area reduction factors due to the ribs between the GC and GDL.
    Jl_agc_agdl_red = Jl.agc_agdl * (pp.Wagc * pp.Lgc) / (pp.Aact / pp.nb_channel_in_gc)
    Jl_cgdl_cgc_red = Jl.cgdl_cgc * (pp.Wcgc * pp.Lgc) / (pp.Aact / pp.nb_channel_in_gc)
    H_gdl_node = pp.Hgdl / NB_GDL   # thickness of one GDL node
    H_mpl_node = pp.Hmpl / NB_MPL   # thickness of one MPL node
    dT(layer,j=0)=temperature_derivative===nothing ? 0.0 :
        (j==0 ? getproperty(temperature_derivative,Symbol(layer,:_T)) : getproperty(temperature_derivative,Symbol(layer,:_T))[j])

    # Differential equations
    #   Anode GDL
    d_s_agdl_dt = ntuple(NB_GDL) do j
        Jl_in = j == 1 ? Jl_agc_agdl_red : Jl.agdl_agdl[j - 1]
        Jl_out = j == NB_GDL ? Jl.agdl_ampl : Jl.agdl_agdl[j]
        liquid_saturation_rate(sv.agdl[j].s,(Jl_in-Jl_out)/H_gdl_node+M_H2O*Sl.agdl[j],pp.epsilon_gdl,T_agdl[j],dT(:agdl,j))
    end

    #   Anode MPL
    d_s_ampl_dt = ntuple(NB_MPL) do j
        Jl_in = j == 1 ? Jl.agdl_ampl : Jl.ampl_ampl[j - 1]
        Jl_out = j == NB_MPL ? Jl.ampl_acl : Jl.ampl_ampl[j]
        liquid_saturation_rate(sv.ampl[j].s,(Jl_in-Jl_out)/H_mpl_node+M_H2O*Sl.ampl[j],pp.epsilon_mpl,T_ampl[j],dT(:ampl,j))
    end

    #   Anode and cathode CLs
    d_s_acl_dt = liquid_saturation_rate(sv.acl.s,
        Jl.ampl_acl/pp.Hacl-M_H2O*S_abs.l_acl+M_H2O*Sl.acl,
        epsilon_cl(:acl,lambda_acl,T_acl,pp.Hacl,pp),T_acl,dT(:acl),cl_porosity_derivative.acl)
    d_s_ccl_dt = liquid_saturation_rate(sv.ccl.s,
        -Jl.ccl_cmpl/pp.Hccl-M_H2O*S_abs.l_ccl+M_H2O*Sl.ccl,
        epsilon_cl(:ccl,lambda_ccl,T_ccl,pp.Hccl,pp),T_ccl,dT(:ccl),cl_porosity_derivative.ccl)

    #   Cathode MPL
    d_s_cmpl_dt = ntuple(NB_MPL) do j
        Jl_in = j == 1 ? Jl.ccl_cmpl : Jl.cmpl_cmpl[j - 1]
        Jl_out = j == NB_MPL ? Jl.cmpl_cgdl : Jl.cmpl_cmpl[j]
        liquid_saturation_rate(sv.cmpl[j].s,(Jl_in-Jl_out)/H_mpl_node+M_H2O*Sl.cmpl[j],pp.epsilon_mpl,T_cmpl[j],dT(:cmpl,j))
    end

    #   Cathode GDL
    d_s_cgdl_dt = ntuple(NB_GDL) do j
        Jl_in = j == 1 ? Jl.cmpl_cgdl : Jl.cgdl_cgdl[j - 1]
        Jl_out = j == NB_GDL ? Jl_cgdl_cgc_red : Jl.cgdl_cgdl[j]
        liquid_saturation_rate(sv.cgdl[j].s,(Jl_in-Jl_out)/H_gdl_node+M_H2O*Sl.cgdl[j],pp.epsilon_gdl,T_cgdl[j],dT(:cgdl,j))
    end

    return MEALiquidWaterDerivative{NB_GDL, NB_MPL}(d_s_agdl_dt, d_s_ampl_dt, d_s_acl_dt, d_s_ccl_dt, d_s_cmpl_dt, d_s_cgdl_dt)
end


"""Calculate the dynamic evolution of water vapour in the porous layers and CLs.

Parameters
----------
sv : CellState1D{NB_GDL, NB_MPL}
    Typed 1D internal state for one gas-channel column.
pp : PhysicalParams
    Fuel-cell physical parameters container (geometry, thicknesses and porous properties).
Jv : MEAVaporFluxes{NB_GDL, NB_MPL}
    Water-vapour inter-layer fluxes (mol·m⁻²·s⁻¹).
Sv : MEAVaporSources{NB_GDL, NB_MPL}
    Vapour phase-change source terms (mol·m⁻³·s⁻¹).
S_abs : MEASorptionSources
    Water absorption/desorption rates at the CL ionomer (mol·m⁻³·s⁻¹).

Returns
-------
CellDerivative1D{NB_GDL, NB_MPL}
    Updated derivative container with C_v derivatives filled in.
"""
function calculate_dyn_vapor_evolution_inside_MEA(
        sv::CellState1D{NB_GDL, NB_MPL},
        pp::PhysicalParams,
        Jv::MEAVaporFluxes{NB_GDL, NB_MPL},
        Sv::MEAVaporSources{NB_GDL, NB_MPL},
        S_abs::MEASorptionSources;
        liquid_derivative::Union{Nothing, MEALiquidWaterDerivative{NB_GDL, NB_MPL}}=nothing,
        cl_porosity_derivative=(acl=0.0, ccl=0.0)
)::MEAVaporDerivative{NB_GDL, NB_MPL} where {NB_GDL, NB_MPL}

    # Extraction of the variables
    s_agc, s_agdl = sv.agc.s, getproperty.(sv.agdl, :s)
    s_ampl, s_acl = getproperty.(sv.ampl, :s), sv.acl.s
    s_ccl, s_cmpl = sv.ccl.s, getproperty.(sv.cmpl, :s)
    s_cgdl, s_cgc = getproperty.(sv.cgdl, :s), sv.cgc.s

    T_agc, T_agdl = sv.agc.T, getproperty.(sv.agdl, :T)
    T_ampl, T_acl = getproperty.(sv.ampl, :T), sv.acl.T
    T_ccl, T_cmpl = sv.ccl.T, getproperty.(sv.cmpl, :T)
    T_cgdl, T_cgc = getproperty.(sv.cgdl, :T), sv.cgc.T

    lambda_acl, lambda_mem, lambda_ccl = sv.acl.lambda, sv.mem.lambda, sv.ccl.lambda
    C_v_agdl, C_v_ampl = getproperty.(sv.agdl, :C_v), getproperty.(sv.ampl, :C_v)
    C_v_cmpl, C_v_cgdl = getproperty.(sv.cmpl, :C_v), getproperty.(sv.cgdl, :C_v)
    ds_agdl = liquid_derivative === nothing ? ntuple(_ -> 0.0, NB_GDL) : liquid_derivative.agdl_s
    ds_ampl = liquid_derivative === nothing ? ntuple(_ -> 0.0, NB_MPL) : liquid_derivative.ampl_s
    ds_acl = liquid_derivative === nothing ? 0.0 : liquid_derivative.acl_s
    ds_ccl = liquid_derivative === nothing ? 0.0 : liquid_derivative.ccl_s
    ds_cmpl = liquid_derivative === nothing ? ntuple(_ -> 0.0, NB_MPL) : liquid_derivative.cmpl_s
    ds_cgdl = liquid_derivative === nothing ? ntuple(_ -> 0.0, NB_GDL) : liquid_derivative.cgdl_s

    # Surface-area reduction factors due to the ribs between the GC and GDL.
    Jv_agc_agdl_red = Jv.agc_agdl * (pp.Wagc * pp.Lgc) / (pp.Aact / pp.nb_channel_in_gc)
    Jv_cgdl_cgc_red = Jv.cgdl_cgc * (pp.Wcgc * pp.Lgc) / (pp.Aact / pp.nb_channel_in_gc)
    H_gdl_node = pp.Hgdl / NB_GDL
    H_mpl_node = pp.Hmpl / NB_MPL

    # Differential equations
    #   Anode GDL
    d_C_v_agdl_dt = ntuple(NB_GDL) do j
        Jv_in = j == 1 ? Jv_agc_agdl_red : Jv.agdl_agdl[j - 1]
        Jv_out = j == NB_GDL ? Jv.agdl_ampl : Jv.agdl_agdl[j]
        gas_concentration_rate(C_v_agdl[j], (Jv_in - Jv_out) / H_gdl_node + Sv.agdl[j],
                               pp.epsilon_gdl, s_agdl[j], ds_agdl[j])
    end

    #   Anode MPL
    d_C_v_ampl_dt = ntuple(NB_MPL) do j
        Jv_in = j == 1 ? Jv.agdl_ampl : Jv.ampl_ampl[j - 1]
        Jv_out = j == NB_MPL ? Jv.ampl_acl : Jv.ampl_ampl[j]
        gas_concentration_rate(C_v_ampl[j], (Jv_in - Jv_out) / H_mpl_node + Sv.ampl[j],
                               pp.epsilon_mpl, s_ampl[j], ds_ampl[j])
    end

    #    Anode and cathode CLs
    d_C_v_acl_dt = gas_concentration_rate(sv.acl.C_v,
        Jv.ampl_acl / pp.Hacl - S_abs.v_acl + Sv.acl,
        epsilon_cl(:acl, lambda_acl, T_acl, pp.Hacl, pp), s_acl, ds_acl,
        cl_porosity_derivative.acl)
    d_C_v_ccl_dt = gas_concentration_rate(sv.ccl.C_v,
        -Jv.ccl_cmpl / pp.Hccl - S_abs.v_ccl + Sv.ccl,
        epsilon_cl(:ccl, lambda_ccl, T_ccl, pp.Hccl, pp), s_ccl, ds_ccl,
        cl_porosity_derivative.ccl)

    #   Cathode MPL
    d_C_v_cmpl_dt = ntuple(NB_MPL) do j
        Jv_in = j == 1 ? Jv.ccl_cmpl : Jv.cmpl_cmpl[j - 1]
        Jv_out = j == NB_MPL ? Jv.cmpl_cgdl : Jv.cmpl_cmpl[j]
        gas_concentration_rate(C_v_cmpl[j], (Jv_in - Jv_out) / H_mpl_node + Sv.cmpl[j],
                               pp.epsilon_mpl, s_cmpl[j], ds_cmpl[j])
    end

    #   Cathode GDL
    d_C_v_cgdl_dt = ntuple(NB_GDL) do j
        Jv_in = j == 1 ? Jv.cmpl_cgdl : Jv.cgdl_cgdl[j - 1]
        Jv_out = j == NB_GDL ? Jv_cgdl_cgc_red : Jv.cgdl_cgdl[j]
        gas_concentration_rate(C_v_cgdl[j], (Jv_in - Jv_out) / H_gdl_node + Sv.cgdl[j],
                               pp.epsilon_gdl, s_cgdl[j], ds_cgdl[j])
    end

    return MEAVaporDerivative{NB_GDL, NB_MPL}(d_C_v_agdl_dt, d_C_v_ampl_dt, d_C_v_acl_dt, d_C_v_ccl_dt, d_C_v_cmpl_dt, d_C_v_cgdl_dt)
end


"""Calculate the dynamic evolution of H₂ (anode) and O₂ (cathode) in the porous layers.

Parameters
----------
sv : CellState1D{NB_GDL, NB_MPL}
    Typed 1D internal state for one gas-channel column.
pp : PhysicalParams
    Fuel-cell physical parameters container (geometry, thicknesses and porous properties).
J_H2 : MEAHydrogenFluxes{NB_GDL, NB_MPL}
    Hydrogen inter-layer fluxes (mol·m⁻²·s⁻¹).
J_O2 : MEAOxygenFluxes{NB_GDL, NB_MPL}
    Oxygen inter-layer fluxes (mol·m⁻²·s⁻¹).
S_H2 : MEAGasReactionSources
    Hydrogen reaction / crossover source terms (mol·m⁻³·s⁻¹).
S_O2 : MEAGasReactionSources
    Oxygen reaction / crossover source terms (mol·m⁻³·s⁻¹).

Returns
-------
CellDerivative1D{NB_GDL, NB_MPL}
    Updated derivative container with C_H2 and C_O2 derivatives filled in.
"""
function calculate_dyn_H2_O2_N2_evolution_inside_MEA(
        sv::CellState1D{NB_GDL, NB_MPL},
        pp::PhysicalParams,
        J_H2::MEAHydrogenFluxes{NB_GDL, NB_MPL},
        J_O2::MEAOxygenFluxes{NB_GDL, NB_MPL},
        J_N2::MEANitrogenFluxes{NB_GDL, NB_MPL},
        S_H2::MEAGasReactionSources,
        S_O2::MEAGasReactionSources;
        multireaction_sources=nothing,
        J_CO2=nothing,
        liquid_derivative::Union{Nothing, MEALiquidWaterDerivative{NB_GDL, NB_MPL}}=nothing,
        cl_porosity_derivative=(acl=0.0, ccl=0.0)
)::MEAGasSpeciesDerivative{NB_GDL, NB_MPL} where {NB_GDL, NB_MPL}

    J_CO2 === nothing && (J_CO2 = MEACarbonDioxideFluxes{NB_GDL,NB_MPL}(0.0,zeros(NB_GDL-1),0.0,zeros(NB_MPL-1),0.0,0.0,zeros(NB_MPL-1),0.0,zeros(NB_GDL-1),0.0))
    # Extraction of the variables
    s_agdl, s_ampl, s_acl = getproperty.(sv.agdl, :s), getproperty.(sv.ampl, :s), sv.acl.s
    s_ccl, s_cmpl, s_cgdl = sv.ccl.s, getproperty.(sv.cmpl, :s), getproperty.(sv.cgdl, :s)

    T_acl, T_ccl = sv.acl.T, sv.ccl.T
    lambda_acl, lambda_ccl = sv.acl.lambda, sv.ccl.lambda
    ds_agdl = liquid_derivative === nothing ? ntuple(_ -> 0.0, NB_GDL) : liquid_derivative.agdl_s
    ds_ampl = liquid_derivative === nothing ? ntuple(_ -> 0.0, NB_MPL) : liquid_derivative.ampl_s
    ds_acl = liquid_derivative === nothing ? 0.0 : liquid_derivative.acl_s
    ds_ccl = liquid_derivative === nothing ? 0.0 : liquid_derivative.ccl_s
    ds_cmpl = liquid_derivative === nothing ? ntuple(_ -> 0.0, NB_MPL) : liquid_derivative.cmpl_s
    ds_cgdl = liquid_derivative === nothing ? ntuple(_ -> 0.0, NB_GDL) : liquid_derivative.cgdl_s

    # Surface-area reduction factors due to the ribs between the GC and GDL.
    J_H2_agc_agdl_red = J_H2.agc_agdl * (pp.Wagc * pp.Lgc) / (pp.Aact / pp.nb_channel_in_gc)
    J_O2_cgdl_cgc_red = J_O2.cgdl_cgc * (pp.Wcgc * pp.Lgc) / (pp.Aact / pp.nb_channel_in_gc)
    J_N2_agc_agdl_red = J_N2.agc_agdl * (pp.Wagc * pp.Lgc) / (pp.Aact / pp.nb_channel_in_gc)
    J_O2_agc_agdl_red = J_O2.agc_agdl * (pp.Wagc * pp.Lgc) / (pp.Aact / pp.nb_channel_in_gc)
    J_N2_cgdl_cgc_red = J_N2.cgdl_cgc * (pp.Wcgc * pp.Lgc) / (pp.Aact / pp.nb_channel_in_gc)
    J_H2_cgdl_cgc_red = J_H2.cgdl_cgc * (pp.Wcgc * pp.Lgc) / (pp.Aact / pp.nb_channel_in_gc)

    H_gdl_node = pp.Hgdl / NB_GDL
    H_mpl_node = pp.Hmpl / NB_MPL

    # Differential equations
    #   Anode GDL
    d_C_H2_agdl_dt = ntuple(NB_GDL) do j
        J_in = j == 1 ? J_H2_agc_agdl_red : J_H2.agdl_agdl[j - 1]
        J_out = j == NB_GDL ? J_H2.agdl_ampl : J_H2.agdl_agdl[j]
        gas_concentration_rate(sv.agdl[j].C_H2, (J_in - J_out) / H_gdl_node,
                               pp.epsilon_gdl, s_agdl[j], ds_agdl[j])
    end
    d_C_N2_agdl_dt = ntuple(NB_GDL) do j
        J_in = j == 1 ? J_N2_agc_agdl_red : J_N2.agdl_agdl[j - 1]
        J_out = j == NB_GDL ? J_N2.agdl_ampl : J_N2.agdl_agdl[j]
        gas_concentration_rate(sv.agdl[j].C_N2, (J_in - J_out) / H_gdl_node,
                               pp.epsilon_gdl, s_agdl[j], ds_agdl[j])
    end
    d_C_O2_agdl_dt = ntuple(NB_GDL) do j
        J_in = j == 1 ? J_O2_agc_agdl_red : J_O2.agdl_agdl[j - 1]
        J_out = j == NB_GDL ? J_O2.agdl_ampl : J_O2.agdl_agdl[j]
        gas_concentration_rate(sv.agdl[j].C_O2, (J_in - J_out) / H_gdl_node,
                               pp.epsilon_gdl, s_agdl[j], ds_agdl[j])
    end

    #   Anode MPL
    d_C_H2_ampl_dt = ntuple(NB_MPL) do j
        J_in = j == 1 ? J_H2.agdl_ampl : J_H2.ampl_ampl[j - 1]
        J_out = j == NB_MPL ? J_H2.ampl_acl : J_H2.ampl_ampl[j]
        gas_concentration_rate(sv.ampl[j].C_H2, (J_in - J_out) / H_mpl_node,
                               pp.epsilon_mpl, s_ampl[j], ds_ampl[j])
    end
    d_C_N2_ampl_dt = ntuple(NB_MPL) do j
        J_in = j == 1 ? J_N2.agdl_ampl : J_N2.ampl_ampl[j - 1]
        J_out = j == NB_MPL ? J_N2.ampl_acl : J_N2.ampl_ampl[j]
        gas_concentration_rate(sv.ampl[j].C_N2, (J_in - J_out) / H_mpl_node,
                               pp.epsilon_mpl, s_ampl[j], ds_ampl[j])
    end
    d_C_O2_ampl_dt = ntuple(NB_MPL) do j
        J_in = j == 1 ? J_O2.agdl_ampl : J_O2.ampl_ampl[j - 1]
        J_out = j == NB_MPL ? J_O2.ampl_acl : J_O2.ampl_ampl[j]
        gas_concentration_rate(sv.ampl[j].C_O2, (J_in - J_out) / H_mpl_node,
                               pp.epsilon_mpl, s_ampl[j], ds_ampl[j])
    end

    #   Anode CL
    S_H2_acl = multireaction_sources === nothing ? -(S_H2.reac + S_H2.cros) : multireaction_sources.anode.H2
    S_O2_acl = multireaction_sources === nothing ? 0.0 : multireaction_sources.anode.O2
    S_H2_ccl = multireaction_sources === nothing ? 0.0 : multireaction_sources.cathode.H2
    S_O2_ccl = multireaction_sources === nothing ? -(S_O2.reac + S_O2.cros) : multireaction_sources.cathode.O2

    epsilon_acl = epsilon_cl(:acl, lambda_acl, T_acl, pp.Hacl, pp)
    d_C_H2_acl_dt = gas_concentration_rate(sv.acl.C_H2,
        J_H2.ampl_acl / pp.Hacl + S_H2_acl, epsilon_acl, s_acl, ds_acl,
        cl_porosity_derivative.acl)
    d_C_N2_acl_dt = gas_concentration_rate(sv.acl.C_N2,
        J_N2.ampl_acl / pp.Hacl, epsilon_acl, s_acl, ds_acl,
        cl_porosity_derivative.acl)
    d_C_O2_acl_dt = gas_concentration_rate(sv.acl.C_O2,
        J_O2.ampl_acl / pp.Hacl + S_O2_acl, epsilon_acl, s_acl, ds_acl,
        cl_porosity_derivative.acl)

    #   Cathode CL
    epsilon_ccl = epsilon_cl(:ccl, lambda_ccl, T_ccl, pp.Hccl, pp)
    d_C_O2_ccl_dt = gas_concentration_rate(sv.ccl.C_O2,
        -J_O2.ccl_cmpl / pp.Hccl + S_O2_ccl, epsilon_ccl, s_ccl, ds_ccl,
        cl_porosity_derivative.ccl)
    d_C_N2_ccl_dt = gas_concentration_rate(sv.ccl.C_N2,
        -J_N2.ccl_cmpl / pp.Hccl, epsilon_ccl, s_ccl, ds_ccl,
        cl_porosity_derivative.ccl)
    d_C_H2_ccl_dt = gas_concentration_rate(sv.ccl.C_H2,
        -J_H2.ccl_cmpl / pp.Hccl + S_H2_ccl, epsilon_ccl, s_ccl, ds_ccl,
        cl_porosity_derivative.ccl)

    #   Cathode MPL
    d_C_O2_cmpl_dt = ntuple(NB_MPL) do j
        J_in = j == 1 ? J_O2.ccl_cmpl : J_O2.cmpl_cmpl[j - 1]
        J_out = j == NB_MPL ? J_O2.cmpl_cgdl : J_O2.cmpl_cmpl[j]
        gas_concentration_rate(sv.cmpl[j].C_O2, (J_in - J_out) / H_mpl_node,
                               pp.epsilon_mpl, s_cmpl[j], ds_cmpl[j])
    end
    d_C_N2_cmpl_dt = ntuple(NB_MPL) do j
        J_in = j == 1 ? J_N2.ccl_cmpl : J_N2.cmpl_cmpl[j - 1]
        J_out = j == NB_MPL ? J_N2.cmpl_cgdl : J_N2.cmpl_cmpl[j]
        gas_concentration_rate(sv.cmpl[j].C_N2, (J_in - J_out) / H_mpl_node,
                               pp.epsilon_mpl, s_cmpl[j], ds_cmpl[j])
    end
    d_C_H2_cmpl_dt = ntuple(NB_MPL) do j
        J_in = j == 1 ? J_H2.ccl_cmpl : J_H2.cmpl_cmpl[j - 1]
        J_out = j == NB_MPL ? J_H2.cmpl_cgdl : J_H2.cmpl_cmpl[j]
        gas_concentration_rate(sv.cmpl[j].C_H2, (J_in - J_out) / H_mpl_node,
                               pp.epsilon_mpl, s_cmpl[j], ds_cmpl[j])
    end

    #   Cathode GDL
    d_C_O2_cgdl_dt = ntuple(NB_GDL) do j
        J_in = j == 1 ? J_O2.cmpl_cgdl : J_O2.cgdl_cgdl[j - 1]
        J_out = j == NB_GDL ? J_O2_cgdl_cgc_red : J_O2.cgdl_cgdl[j]
        gas_concentration_rate(sv.cgdl[j].C_O2, (J_in - J_out) / H_gdl_node,
                               pp.epsilon_gdl, s_cgdl[j], ds_cgdl[j])
    end
    d_C_N2_cgdl_dt = ntuple(NB_GDL) do j
        J_in = j == 1 ? J_N2.cmpl_cgdl : J_N2.cgdl_cgdl[j - 1]
        J_out = j == NB_GDL ? J_N2_cgdl_cgc_red : J_N2.cgdl_cgdl[j]
        gas_concentration_rate(sv.cgdl[j].C_N2, (J_in - J_out) / H_gdl_node,
                               pp.epsilon_gdl, s_cgdl[j], ds_cgdl[j])
    end
    d_C_H2_cgdl_dt = ntuple(NB_GDL) do j
        J_in = j == 1 ? J_H2.cmpl_cgdl : J_H2.cgdl_cgdl[j - 1]
        J_out = j == NB_GDL ? J_H2_cgdl_cgc_red : J_H2.cgdl_cgdl[j]
        gas_concentration_rate(sv.cgdl[j].C_H2, (J_in - J_out) / H_gdl_node,
                               pp.epsilon_gdl, s_cgdl[j], ds_cgdl[j])
    end

    J_CO2_agc_agdl_red = J_CO2.agc_agdl * (pp.Wagc * pp.Lgc) / (pp.Aact / pp.nb_channel_in_gc)
    J_CO2_cgdl_cgc_red = J_CO2.cgdl_cgc * (pp.Wcgc * pp.Lgc) / (pp.Aact / pp.nb_channel_in_gc)
    d_C_CO2_agdl_dt = ntuple(NB_GDL) do j
        J_in = j == 1 ? J_CO2_agc_agdl_red : J_CO2.agdl_agdl[j - 1]
        J_out = j == NB_GDL ? J_CO2.agdl_ampl : J_CO2.agdl_agdl[j]
        gas_concentration_rate(sv.agdl[j].C_CO2, (J_in - J_out) / H_gdl_node,
                               pp.epsilon_gdl, s_agdl[j], ds_agdl[j])
    end
    d_C_CO2_ampl_dt = ntuple(NB_MPL) do j
        J_in = j == 1 ? J_CO2.agdl_ampl : J_CO2.ampl_ampl[j - 1]
        J_out = j == NB_MPL ? J_CO2.ampl_acl : J_CO2.ampl_ampl[j]
        gas_concentration_rate(sv.ampl[j].C_CO2, (J_in - J_out) / H_mpl_node,
                               pp.epsilon_mpl, s_ampl[j], ds_ampl[j])
    end
    d_C_CO2_acl_dt = gas_concentration_rate(sv.acl.C_CO2,
        J_CO2.ampl_acl / pp.Hacl + (multireaction_sources === nothing ? 0.0 : multireaction_sources.anode.CO2), epsilon_acl, s_acl, ds_acl,
        cl_porosity_derivative.acl)
    d_C_CO2_ccl_dt = gas_concentration_rate(sv.ccl.C_CO2,
        -J_CO2.ccl_cmpl / pp.Hccl + (multireaction_sources === nothing ? 0.0 : multireaction_sources.cathode.CO2), epsilon_ccl, s_ccl, ds_ccl,
        cl_porosity_derivative.ccl)
    d_C_CO2_cmpl_dt = ntuple(NB_MPL) do j
        J_in = j == 1 ? J_CO2.ccl_cmpl : J_CO2.cmpl_cmpl[j - 1]
        J_out = j == NB_MPL ? J_CO2.cmpl_cgdl : J_CO2.cmpl_cmpl[j]
        gas_concentration_rate(sv.cmpl[j].C_CO2, (J_in - J_out) / H_mpl_node,
                               pp.epsilon_mpl, s_cmpl[j], ds_cmpl[j])
    end
    d_C_CO2_cgdl_dt = ntuple(NB_GDL) do j
        J_in = j == 1 ? J_CO2.cmpl_cgdl : J_CO2.cgdl_cgdl[j - 1]
        J_out = j == NB_GDL ? J_CO2_cgdl_cgc_red : J_CO2.cgdl_cgdl[j]
        gas_concentration_rate(sv.cgdl[j].C_CO2, (J_in - J_out) / H_gdl_node,
                               pp.epsilon_gdl, s_cgdl[j], ds_cgdl[j])
    end
    return MEAGasSpeciesDerivative{NB_GDL, NB_MPL}(
        d_C_H2_agdl_dt, d_C_H2_ampl_dt, d_C_H2_acl_dt,
        d_C_N2_agdl_dt, d_C_N2_ampl_dt, d_C_N2_acl_dt,
        d_C_O2_ccl_dt, d_C_O2_cmpl_dt, d_C_O2_cgdl_dt,
        d_C_N2_ccl_dt, d_C_N2_cmpl_dt, d_C_N2_cgdl_dt,
        d_C_O2_agdl_dt, d_C_O2_ampl_dt, d_C_O2_acl_dt,
        d_C_H2_ccl_dt, d_C_H2_cmpl_dt, d_C_H2_cgdl_dt,
        d_C_CO2_agdl_dt,d_C_CO2_ampl_dt,d_C_CO2_acl_dt,d_C_CO2_ccl_dt,d_C_CO2_cmpl_dt,d_C_CO2_cgdl_dt
    )
end

"""Evaluate one column's multireaction kinetics, gas/water sources and potentials."""
function calculate_multireaction_coupling(sv::CellState1D, i_fc::Real,
                                          C_O2_Pt::Real, fc::AbstractFuelCell,
                                          cfg::SimulationConfig)
    pp = fc.physical_parameters
    rp = cfg.reaction_parameters
    coverage(cl) = cfg.enable_pt_oxide ? PtCoverage(cl.theta_PtOH, cl.theta_Pt_sO, cl.theta_Pt_bO) : nothing
    anode = ElectrodeState(sv.acl.T, sv.acl.C_H2, sv.acl.C_O2, sv.acl.phi_a, pp.Hacl, coverage(sv.acl),sv.acl.C_CO2,cfg.enable_cor)
    cathode = ElectrodeState(sv.ccl.T, sv.ccl.C_H2, sv.ccl.C_O2, sv.ccl.phi_c, pp.Hccl, coverage(sv.ccl),sv.ccl.C_CO2,cfg.enable_cor)
    reactions = multireaction(anode, cathode, rp; C_O2_Pt=C_O2_Pt)
    potentials = potential_derivatives(reactions, i_fc, rp)

    # Signed membrane fluxes are positive from ACL to CCL. Retaining both
    # boundary concentrations lets crossover vanish at equilibrium and reverse
    # when the concentration gradient reverses.
    J_H2_mem = k_H2(sv.mem.lambda, sv.mem.T, pp.kappa_co, pp) * R * sv.mem.T /
               pp.Hmem * (sv.acl.C_H2 - sv.ccl.C_H2)
    J_O2_mem = k_O2(sv.mem.lambda, sv.mem.T, pp.kappa_co, pp) * R * sv.mem.T /
               pp.Hmem * (sv.acl.C_O2 - sv.ccl.C_O2)
    sources = reaction_sources(reactions, pp.Hacl, pp.Hccl;
                               J_H2=J_H2_mem, J_O2=J_O2_mem)
    water = MEAWaterProductionSources(sources.anode.H2O, sources.cathode.H2O)
    anode_heat = reaction_heat(anode, reactions.anode, rp)
    cathode_heat = reaction_heat(cathode, reactions.cathode, rp)
    heat = MEAReactionHeat(anode_heat.total, cathode_heat.total)
    voltage = MEAVoltageDerivative(0.0, potentials.phi_a, potentials.phi_c)
    rates(e) = e.oxide === nothing ? (0.0, 0.0, 0.0) :
        (e.oxide.rates.OH, e.oxide.rates.sO, e.oxide.rates.bO)
    pt_rates = (rates(reactions.anode), rates(reactions.cathode))
    return (; reactions, sources, water, heat, voltage, pt_rates)
end


"""Calculate the dynamic evolution of the cathode overpotential eta_c.

Parameters
----------
i_fc
    Fuel cell current density (A·m⁻²).
C_O2_Pt
    Oxygen concentration at the platinum surface (mol·m⁻³).
T_ccl : Float64
    Temperature in the cathode catalyst layer (K).
eta_c : Float64
    Cathode overpotential (V).
pp : PhysicalParams
    Fuel-cell physical parameters container.
i_n
    Crossover current density (A·m⁻²).

Returns
-------
CellDerivative1D{NB_GDL, NB_MPL}
    Updated derivative container with eta_c derivative filled in (ccl).
"""
function calculate_dyn_voltage_evolution(
        i_fc,
        C_O2_Pt,
        T_ccl::Float64,
        eta_c::Float64,
        pp::PhysicalParams,
        i_n
)::MEAVoltageDerivative

    # Extraction of the parameters
    alpha_c = pp.alpha_c # Charge transfer coefficient for the cathode reaction

    # During nonlinear/DAE iterations the algebraic unknown `C_O2_Pt` can briefly
    # step outside its physical domain (e.g. slightly negative). Protect the
    # Butler–Volmer concentration term which uses fractional powers.
    C_O2_Pt_safe = _positive_concentration_value(C_O2_Pt)

    # Differential equation
    d_eta_c_ccl_dt = 1 / (pp.C_scl * pp.Hccl) * ((i_fc + i_n) -
             pp.i0_c_ref * (C_O2_Pt_safe / C_O2ref_red)^pp.kappa_c *
             exp(-Eact_O2_red / R * (1 / T_ccl - 1 / Tref_O2_red)) *
             exp(alpha_c * F / (R * T_ccl) * eta_c))

    return MEAVoltageDerivative(d_eta_c_ccl_dt)
end


"""Calculate the dynamic evolution of temperature throughout the MEA.

Parameters
----------
pp : PhysicalParams
    Fuel-cell physical parameters container (layer thicknesses come from `pp`).
nb_gdl, nb_mpl : Int64
    Node counts (match NB_GDL and NB_MPL).
rho_Cp0 : MEAThermalIntermediates{NB_GDL, NB_MPL}
    Volumetric heat capacities at each node (J·m⁻³·K⁻¹).
Jt : MEAThermalFluxes{NB_GDL, NB_MPL}
    Conductive thermal fluxes through the MEA (W·m⁻²).
Q_r : MEAReactionHeat
    Electrochemical reaction heat sources (W·m⁻³).
Q_sorp : MEASorptionHeat
    Sorption heat sources (W·m⁻³).
Q_liq : MEALiquidHeat{NB_GDL, NB_MPL}
    Liquefaction / evaporation heat sources (W·m⁻³).
Q_p : MEAProtonHeat
    Ionic (protonic) Joule-heating sources (W·m⁻³).
Q_e : MEAElectricHeat{NB_GDL, NB_MPL}
    Electronic Joule-heating sources (W·m⁻³).

Returns
-------
CellDerivative1D{NB_GDL, NB_MPL}
    Updated derivative container with T derivatives filled in for all MEA layers.
"""
function calculate_dyn_temperature_evolution_inside_MEA(
        rho_Cp0::MEAThermalIntermediates{NB_GDL, NB_MPL},
        pp::PhysicalParams,
        Jt::MEAThermalFluxes{NB_GDL, NB_MPL},
        Q_r::MEAReactionHeat,
        Q_sorp::MEASorptionHeat,
        Q_liq::MEALiquidHeat{NB_GDL, NB_MPL},
        Q_p::MEAProtonHeat,
        Q_e::MEAElectricHeat{NB_GDL, NB_MPL}
)::MEATemperatureDerivative{NB_GDL, NB_MPL} where {NB_GDL, NB_MPL}

    # Extraction of the parameters
    H_gdl_node = pp.Hgdl / NB_GDL
    H_mpl_node = pp.Hmpl / NB_MPL

    # Differential equations
    #   Anode GDL
    d_T_agdl_dt = ntuple(NB_GDL) do j
        Jt_in = j == 1 ? Jt.agc_agdl : Jt.agdl_agdl[j - 1]
        Jt_out = j == NB_GDL ? Jt.agdl_ampl : Jt.agdl_agdl[j]
        (1 / rho_Cp0.agdl[j]) * ((Jt_in - Jt_out) / H_gdl_node + Q_liq.agdl[j] + Q_e.agdl[j])
    end

    #   Anode MPL
    d_T_ampl_dt = ntuple(NB_MPL) do j
        Jt_in = j == 1 ? Jt.agdl_ampl : Jt.ampl_ampl[j - 1]
        Jt_out = j == NB_MPL ? Jt.ampl_acl : Jt.ampl_ampl[j]
        (1 / rho_Cp0.ampl[j]) * ((Jt_in - Jt_out) / H_mpl_node + Q_liq.ampl[j] + Q_e.ampl[j])
    end

    #  Anode and cathode CLs + membrane
    d_T_acl_dt = (1 / rho_Cp0.acl) * ((Jt.ampl_acl - Jt.acl_mem) / pp.Hacl +
            Q_r.acl + Q_sorp.v_acl + Q_sorp.l_acl + Q_liq.acl + Q_e.acl)
    d_T_mem_dt = (1 / rho_Cp0.mem) * ((Jt.acl_mem - Jt.mem_ccl) / pp.Hmem + Q_p.mem)
    d_T_ccl_dt = (1 / rho_Cp0.ccl) * ((Jt.mem_ccl - Jt.ccl_cmpl) / pp.Hccl +
            Q_r.ccl + Q_sorp.v_ccl + Q_sorp.l_ccl + Q_liq.ccl + Q_p.ccl + Q_e.ccl)

    #   Cathode MPL
    d_T_cmpl_dt = ntuple(NB_MPL) do j
        Jt_in = j == 1 ? Jt.ccl_cmpl : Jt.cmpl_cmpl[j - 1]
        Jt_out = j == NB_MPL ? Jt.cmpl_cgdl : Jt.cmpl_cmpl[j]
        (1 / rho_Cp0.cmpl[j]) * ((Jt_in - Jt_out) / H_mpl_node + Q_liq.cmpl[j] + Q_e.cmpl[j])
    end

    #  Cathode GDL
    d_T_cgdl_dt = ntuple(NB_GDL) do j
        Jt_in = j == 1 ? Jt.cmpl_cgdl : Jt.cgdl_cgdl[j - 1]
        Jt_out = j == NB_GDL ? Jt.cgdl_cgc : Jt.cgdl_cgdl[j]
        (1 / rho_Cp0.cgdl[j]) * ((Jt_in - Jt_out) / H_gdl_node + Q_liq.cgdl[j] + Q_e.cgdl[j])
    end

    return MEATemperatureDerivative{NB_GDL, NB_MPL}(d_T_agdl_dt, d_T_ampl_dt, d_T_acl_dt, d_T_mem_dt, d_T_ccl_dt, d_T_cmpl_dt, d_T_cgdl_dt)
end
