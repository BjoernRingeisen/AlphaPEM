using Test
using AlphaPEM

@testset "Reaction heat thermochemical basis" begin
    C, E, U = AlphaPEM.Config, AlphaPEM.Core.Electrochemistry, AlphaPEM.Utils
    T, F = 298.15, U.F
    p = C.ReactionParams()
    c_ref = U.Pref_eq / (U.R * T)
    # Prescribe charge-balanced half-reaction currents to isolate heat from
    # kinetic fitting, inventory limiting, Pt oxidation and corrosion.
    reaction(j, equilibrium, eta) = E.ReactionResult(j, j, equilibrium, eta, 1.0, 1.0, c_ref)
    for j in (1e4, -1e4), (Ha, Hc) in ((3e-6, 13e-6), (1e-5, 2e-5)),
        (eta_h, eta_o) in ((0.0, 0.0), (0.05, 0.3))
        a = E.ElectrodeState(T, c_ref, c_ref, eta_h, Ha)
        c = E.ElectrodeState(T, c_ref, c_ref, U.E0 - eta_o, Hc)
        ra = E.ElectrodeReactions(reaction(j, 0.0, eta_h), reaction(0, U.E0, 0.0), j)
        rc = E.ElectrodeReactions(reaction(0, 0.0, 0.0), reaction(j, U.E0, eta_o), -j)
        ha, hc = E.reaction_heat(a, ra, p), E.reaction_heat(c, rc, p)
        q_area = Ha * ha.total + Hc * hc.total
        voltage = c.phi - a.phi
        thermoneutral = voltage + q_area / j
        @test ha.nernst == 0.0 && hc.nernst == 0.0
        @test thermoneutral ≈ 1.48114655 atol=5e-9
        # Independent HHV benchmark: NIST liquid-water formation enthalpy.
        # Rounded source entropies/E0 leave a 0.063 mV residual, not exact equality.
        @test thermoneutral ≈ 285830.0 / (2F) atol=7e-5
        sources = E.reaction_sources(E.MultireactionResult(ra, rc), Ha, Hc)
        @test Ha * sources.anode.H2 ≈ -j / (2F)
        @test Hc * sources.cathode.O2 ≈ -j / (4F)
        @test Hc * sources.cathode.H2O ≈ j / (2F)
    end

    @testset "Legacy heat uses the same two-electron entropy extent" begin
        M, Types = AlphaPEM.Core.Models, AlphaPEM.Core.Types
        np = C.NumericalParams(nb_gc=1, nb_gdl=2, nb_mpl=1)
        cfg = C.SimulationConfig(type_display=:no_display, type_auxiliary=:no_auxiliary,
                                 numerical_parameters=np)
        fc = AlphaPEM.Fuelcell.create_fuelcell(:ZSW_nominal, :full; year=2024, nb_gc=1)
        sim = M.AlphaPEM(fc, AlphaPEM.Currents.create_current(cfg.type_current, fc), cfg)
        names = M.canonical_cell_solver_variable_names_1D(np.nb_gdl, np.nb_mpl)
        initial = M.create_initial_variable_values(sim)
        state = M._unpack_cell_state_1D(initial[1:length(names)], np.nb_gdl, np.nb_mpl)
        j = 1e4
        flows = M.calculate_flows_1D_MEA!(Types.MEAFlowsWorkspace(np.nb_gdl, np.nb_mpl),
            Types.MEAFlowsIntWorkspace(np.nb_gdl, np.nb_mpl), state, j, 1.0, 1.0, fc, cfg)
        heat = M.calculate_heat_transfers!(Types.MEAHeatWorkspace(np.nb_gdl, np.nb_mpl),
            Types.MEAHeatIntWorkspace(np.nb_gdl, np.nb_mpl), state, j, fc, cfg, flows.S_abs, flows.Sl)
        pp = fc.physical_parameters
        @test pp.Hacl * heat.Q_r.acl ≈ -j * state.acl.T * U.delta_s_HOR / (2F)
        @test pp.Hccl * heat.Q_r.ccl - j * state.ccl.eta_c ≈ -j * state.ccl.T * U.delta_s_ORR / (2F)
    end
end

@testset "COR heat follows the reduction-potential convention" begin
    C, E, U = AlphaPEM.Config, AlphaPEM.Core.Electrochemistry, AlphaPEM.Utils

    # Prescribed COR current isolates thermodynamics from exchange-current
    # prefactors, Pt oxide kinetics and the HOR/ORR branches.
    function isolated_cor_heat(temp, activity, j, eta, H, p, catalyzed)
        co2 = activity * U.Pref_eq / (U.R * temp)
        eq = E.reversible_cor_potential(co2, temp, p; catalyzed)
        state = E.ElectrodeState(temp, 0.0, 0.0, eq + eta, H,
                                E.PtCoverage(), co2, true)
        gas_off = E.ReactionResult(0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0)
        active = E.ReactionResult(j, j, eq, eta, 1.0, 1.0, co2)
        inactive = E.ReactionResult(0.0, 0.0, eq, 0.0, 1.0, 1.0, co2)
        n = catalyzed ? 3 : 4
        cor = E.CarbonOxidationResult(catalyzed ? inactive : active,
            catalyzed ? active : inactive, p.carbon_molar_mass * j / (n * U.F), 0.0)
        reactions = E.ElectrodeReactions(gas_off, gas_off, j, nothing, cor)
        return E.reaction_heat(state, reactions, p), state
    end

    p = C.ReactionParams()
    h, _ = isolated_cor_heat(353.15, 1.0, 1000.0, 0.0, 1e-5, p, false)
    coefficient = h.entropy * 1e-5 / 1000.0
    @test coefficient ≈ -0.30128457 atol=5e-9
    @test 2abs(coefficient) ≈ 0.60256915 atol=5e-9

    # Nonzero test entropies ensure the catalyzed branch cannot pass merely
    # because its production entropy coefficient is zero.
    for catalyzed in (false, true), ds_cat in (-125.0, 0.0, 125.0),
        temp in (298.15, 353.15), activity in (0.1, 1.0, 3.0),
        H in (3e-6, 2e-5), eta in (0.0, 0.2)
        p = C.ReactionParams(delta_s_cor_cat=ds_cat)
        j = 1000.0
        h, state = isolated_cor_heat(temp, activity, j, eta, H, p, catalyzed)
        n = catalyzed ? 3 : 4
        dh = catalyzed ? p.delta_h_cor_cat : p.delta_h_cor_noncat
        # Independent temperature derivative of the actual reduction-potential
        # implementation, with gas ACTIVITY held fixed, not concentration.
        potential(T) = E.reversible_cor_potential(
            activity * U.Pref_eq / (U.R * T), T, p; catalyzed)
        dt = 0.01
        slope = (potential(temp + dt) - potential(temp - dt)) / (2dt)
        @test H * (h.entropy + h.nernst) / j ≈ temp * slope atol=1e-9 rtol=1e-8
        @test H * h.activation / j ≈ eta atol=1e-14
        # With constant reduction enthalpy and reference offset, entropy and
        # activity must cancel from total heat at fixed electrode potential.
        @test H * h.total / j ≈ state.phi + p.SHE_absolute_potential + dh / (n * U.F) atol=1e-12
        off, _ = isolated_cor_heat(temp, activity, 0.0, eta, H, p, catalyzed)
        @test off.total == 0.0
    end

    # Directional heat law only: production COR kinetics remain irreversible.
    # Reversing both current and overpotential reverses reversible heat while
    # keeping activation dissipation positive in either direction.
    for catalyzed in (false, true)
        p = C.ReactionParams(delta_s_cor_cat=-125.0)
        forward, _ = isolated_cor_heat(353.15, 2.0, 1000.0, 0.2, 1e-5, p, catalyzed)
        reverse, _ = isolated_cor_heat(353.15, 2.0, -1000.0, -0.2, 1e-5, p, catalyzed)
        @test reverse.entropy ≈ -forward.entropy
        @test reverse.nernst ≈ -forward.nernst
        @test reverse.activation ≈ forward.activation && reverse.activation > 0
    end
end
