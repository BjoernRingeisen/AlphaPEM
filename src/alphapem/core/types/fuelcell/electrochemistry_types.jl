"""Surface PtOH/PtO fractions and bulk oxide in equivalent layers.

Free surface Pt is exactly `1 - OH - sO`; it is not an independent state.
"""
Base.@kwdef struct PtCoverage
    OH::Float64 = 0.0
    sO::Float64 = 0.0
    bO::Float64 = 0.0
end

"""Geometric oxidation currents, conservative coverage rates and water sink."""
struct PtOxideResult
    currents::NTuple{3,Float64}
    rates::PtCoverage
    surface_inventory::Float64
    max_layers::Float64
end

"""Local reaction inputs, with potential vs SHE and optional oxide coverage."""
Base.@kwdef struct ElectrodeState
    T::Float64
    C_H2::Float64
    C_O2::Float64
    phi::Float64
    thickness::Float64
    coverage::Union{Nothing,PtCoverage} = nothing
    C_CO2::Float64 = 0.0
    enable_cor::Bool = false
end

ElectrodeState(T, H2, O2, phi, H) = ElectrodeState(T, H2, O2, phi, H, nothing)
ElectrodeState(T, H2, O2, phi, H, coverage) =
    ElectrodeState(T, H2, O2, phi, H, coverage, 0.0, false)

"""One geometric-area reaction current and its diagnostic factors.

`raw_current` includes concentration availability and precedes the inventory
cap. `concentration` records the gas concentration used for the Nernst potential
and kinetics, so heat uses the same reaction state (especially cathode Pt O₂).
"""
struct ReactionResult
    current::Float64
    raw_current::Float64
    reversible_potential::Float64
    overpotential::Float64
    availability::Float64
    inventory_factor::Float64
    concentration::Float64
end

"""Non-catalyzed and PtOH-catalyzed carbon oxidation (geometric currents)."""
struct CarbonOxidationResult
    noncat::ReactionResult
    cat::ReactionResult
    carbon_loss_rate::Float64 # kg carbon / m² geometric area / s
    ptOH_consumption_rate::Float64 # surface fraction / s
end

"""Local gas, carbon and oxide results with oxidation-positive net Faraday current."""
struct ElectrodeReactions
    hor::ReactionResult
    orr::ReactionResult
    faraday::Float64
    oxide::Union{Nothing,PtOxideResult}
    cor::Union{Nothing,CarbonOxidationResult}
end

ElectrodeReactions(hor, orr, faraday) = ElectrodeReactions(hor, orr, faraday, nothing)
ElectrodeReactions(hor, orr, faraday, oxide) = ElectrodeReactions(hor, orr, faraday, oxide, nothing)

"""Local reaction results at both electrodes of one cell column."""
struct MultireactionResult
    anode::ElectrodeReactions
    cathode::ElectrodeReactions
end
