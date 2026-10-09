# -*- coding: utf-8 -*-

"""This file represents the equations for calculating the cell voltage. It is a component of the fuel cell model.
"""

# _____________________________________________________Cell voltage_____________________________________________________

"""Calculate the cell voltage in volt.

Parameters
----------
i_fc : Float64
    The current density (A/m²).
C_O2_Pt : Float64
    The oxygen concentration at the platinum surface in the cathode catalyst layer (mol/m³).
sv : CellState1D
    The typed 1D cell-column state (MEA+GC) for one gas-channel position.
fc : AbstractFuelCell
    The fuel cell instance providing model parameters.

Returns
-------
Float64
    The cell voltage in volt.
"""
function calculate_cell_voltage(i_fc::Real, C_O2_Pt::Real, sv::CellState1D, fc::AbstractFuelCell)

    # Extraction of the variables
    T_ampl, T_acl = getproperty.(sv.ampl, :T), _positive_temperature_value(sv.acl.T)
    T_mem = _positive_temperature_value(sv.mem.T)
    T_ccl, T_cmpl = _positive_temperature_value(sv.ccl.T), getproperty.(sv.cmpl, :T)

    lambda_mem, lambda_ccl = sv.mem.lambda, sv.ccl.lambda

    C_H2_ampl, C_H2_acl = getproperty.(sv.ampl, :C_H2), _nonnegative_value(sv.acl.C_H2)
    C_O2_ccl, C_O2_cmpl = sv.ccl.C_O2, getproperty.(sv.cmpl, :C_O2)

    eta_c = sv.ccl.eta_c
    C_O2_Pt_safe = _positive_concentration_value(C_O2_Pt)

    # Extraction of the parameters
    pp = fc.physical_parameters
    Hmem, Hacl, Hccl = pp.Hmem, pp.Hacl, pp.Hccl
    Re, kappa_co = pp.Re, pp.kappa_co

    # The equilibrium potential
    Ueq = E0 - 8.5e-4 * (T_ccl - 298.15) + R * T_ccl / (2 * F) *
          (log(R * T_acl * C_H2_acl / Pref_eq) +
           0.5 * log(R * T_ccl * C_O2_Pt_safe / Pref_eq))

    # The crossover current density
    T_acl_mem_ccl = average([T_acl, T_mem, T_ccl],
                            [Hacl / (Hacl + Hmem + Hccl), Hmem / (Hacl + Hmem + Hccl), Hccl / (Hacl + Hmem + Hccl)])
    i_H2 = 2 * F * R * T_acl_mem_ccl / Hmem * C_H2_acl * k_H2(lambda_mem, T_mem, kappa_co, pp)
    i_O2 = 4 * F * R * T_acl_mem_ccl / Hmem * C_O2_ccl * k_O2(lambda_mem, T_mem, kappa_co, pp)
    i_n = i_H2 + i_O2

    # The proton resistance
    #       The proton resistance at the membrane : Rmem
    Rmem = Hmem / sigma_p_eff(:mem, lambda_mem, T_mem, nothing, pp)
    #       The proton resistance at the cathode catalyst layer : Rccl
    Rccl = Hccl / sigma_p_eff(:ccl, lambda_ccl, T_ccl, Hccl, pp)
    #       The total proton resistance
    Rp = Rmem + Rccl  # Its value is around [4-7]e-6 ohm.m².

    # The cell voltage
    Ucell = Ueq - eta_c - (i_fc + i_n) * (Rp + Re)

    return Ucell
end

"""Cell voltage obtained from the two local electrode-potential states.

The potential difference supplies the electrochemical voltage. The external
terminal voltage additionally includes the same membrane, CCL protonic and
electronic losses used by the legacy closure.
"""
function calculate_cell_voltage_from_potentials(i_fc::Real, sv::CellState1D,
                                                fc::AbstractFuelCell)
    pp = fc.physical_parameters
    T_acl = _positive_temperature_value(sv.acl.T)
    T_mem = _positive_temperature_value(sv.mem.T)
    T_ccl = _positive_temperature_value(sv.ccl.T)
    lambda_mem, lambda_ccl = sv.mem.lambda, sv.ccl.lambda
    C_H2_acl = _nonnegative_value(sv.acl.C_H2)
    C_H2_ccl = _nonnegative_value(sv.ccl.C_H2)
    C_O2_acl = _nonnegative_value(sv.acl.C_O2)
    C_O2_ccl = _nonnegative_value(sv.ccl.C_O2)

    T_mean = average([T_acl, T_mem, T_ccl],
                     [pp.Hacl, pp.Hmem, pp.Hccl] ./ (pp.Hacl + pp.Hmem + pp.Hccl))
    # Match the signed membrane-gradient crossover used by the multireaction
    # species balances. Equal inventories on both sides (air/air startup) have
    # zero crossover current instead of the legacy one-way cathode-O2 loss.
    i_H2 = 2F * R * T_mean / pp.Hmem * (C_H2_acl - C_H2_ccl) *
           k_H2(lambda_mem, T_mem, pp.kappa_co, pp)
    i_O2 = 4F * R * T_mean / pp.Hmem * (C_O2_ccl - C_O2_acl) *
           k_O2(lambda_mem, T_mem, pp.kappa_co, pp)
    i_n = i_H2 + i_O2

    Rmem = pp.Hmem / sigma_p_eff(:mem, lambda_mem, T_mem, nothing, pp)
    Rccl = pp.Hccl / (3 * sigma_p_eff(:ccl, lambda_ccl, T_ccl, pp.Hccl, pp))
    return sv.ccl.phi_c - sv.acl.phi_a - (i_fc + i_n) * (Rmem + Rccl + pp.Re)
end
