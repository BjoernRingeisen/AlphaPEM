# -*- coding: utf-8 -*-

"""This file represents the equations for calculating the current distribution through the GC.
It is a component of the fuel cell model.
"""

# _________________________________________________Current distribution_________________________________________________

"""Fill algebraic residuals for the GC current-distribution block in physical units.

Unknown layout is fixed as:
    [U_cell, i_fc[1:nb_gc], C_O2_Pt[1:nb_gc]]
"""
function gc_current_distribution_residuals!(res::AbstractVector{Float64},
                                            U_cell::Float64,
                                            i_fc::AbstractVector{Float64},
                                            C_O2_Pt::AbstractVector{Float64},
                                            i_fc_cell::Float64,
                                            sv::AbstractVector{<:CellState1D},
                                            fc::AbstractFuelCell,
                                            cfg::SimulationConfig)
    nb_gc = cfg.numerical_parameters.nb_gc
    length(sv) == nb_gc || throw(ArgumentError("sv size mismatch with cfg.numerical_parameters.nb_gc in gc_current_distribution_residuals!."))
    length(i_fc) == nb_gc || throw(ArgumentError("i_fc size mismatch in gc_current_distribution_residuals!."))
    length(C_O2_Pt) == nb_gc || throw(ArgumentError("C_O2_Pt size mismatch in gc_current_distribution_residuals!."))
    length(res) == 2 * nb_gc + 1 || throw(ArgumentError("res size mismatch in gc_current_distribution_residuals!."))

    @inbounds for i in 1:nb_gc
        voltage = if cfg.electrochemistry_model == :multireaction_potential
            calculate_cell_voltage_from_potentials(i_fc[i], sv[i], fc)
        else
            calculate_cell_voltage(i_fc[i], C_O2_Pt[i], sv[i], fc)
        end
        res[i] = voltage - U_cell
    end
    res[nb_gc + 1] = i_fc_cell - average(i_fc)
    @inbounds for i in 1:nb_gc
        oxygen_consuming_current = cfg.electrochemistry_model == :multireaction_potential ?
                                   max(i_fc[i], 0.0) : i_fc[i]
        if cfg.enable_pt_oxide
            cl = sv[i].ccl
            state = ElectrodeState(cl.T, cl.C_H2, cl.C_O2, cl.phi_c, fc.physical_parameters.Hccl,
                PtCoverage(cl.theta_PtOH, cl.theta_Pt_sO, cl.theta_Pt_bO),cl.C_CO2,cfg.enable_cor)
            rp = cfg.reaction_parameters
            reaction = electrode_reactions(state, rp.roughness_c, rp; C_O2_Pt=C_O2_Pt[i])
            oxygen_consuming_current = max(reaction.orr.current, 0.0)
        end
        res[nb_gc + 1 + i] = calculate_C_O2_Pt(oxygen_consuming_current, sv[i], fc) - C_O2_Pt[i]
    end
    return nothing
end

