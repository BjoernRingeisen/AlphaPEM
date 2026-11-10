# cell_state.jl
#
# Defines the Julia structs representing the internal state at each spatial node
# of the 1D MEA model and the 1D+1D (P2D) gas channel discretisation.
#
# Model layout (anode → cathode, along the through-plane direction):
#
#   AGC | AGDL×nb_gdl | AMPL×nb_mpl | ACL | Membrane | CCL | CMPL×nb_mpl | CGDL×nb_gdl | CGC
#
# Each struct holds only the state variables that are physically meaningful at
# that layer. In particular:
#   - `lambda` (ionomer water content) only appears in ACL, Membrane and CCL.
#   - `s`      (liquid water saturation) appears in all porous layers and GCs,
#               but the dominant transport mechanism differs:
#               convective in the GC, capillary-diffusive in GDL/MPL/CL.
#   - `C_N2`   is carried by both GC states (always present in AnodeGCState for
#               future flow-through-anode mode; always present in CathodeGCState
#               for air operation).
#   - `eta_c`  (cathode overpotential) is localised at the CCL, where the
#               oxygen reduction reaction takes place.
#   - `phi_a` and `phi_c` are the electrode potentials versus SHE, localised at
#               the ACL and CCL for the multireaction double-layer balances.

# ────────────────────────────────────────────────────────────────────────────────
# Abstract root type
# ────────────────────────────────────────────────────────────────────────────────

abstract type AbstractCellState end

# ────────────────────────────────────────────────────────────────────────────────
# Anode side
# ────────────────────────────────────────────────────────────────────────────────

"""Internal state at one anode gas-channel location."""
struct AnodeGCState <: AbstractCellState
    T    :: Float64   # Temperature                       (K)
    C_v  :: Float64   # Water vapour concentration        (mol·m⁻³)
    s    :: Float64   # Liquid water saturation           (–)
    C_H2 :: Float64   # Hydrogen concentration            (mol·m⁻³)
    C_N2 :: Float64   # Nitrogen concentration            (mol·m⁻³)
    C_O2::Float64   # Additional gas inventory (mol/m^3)
    C_CO2::Float64 # CO₂ gas inventory
end

"""Internal state at one anode gas-diffusion-layer location."""
struct AnodeGDLState <: AbstractCellState
    T    :: Float64   # Temperature                       (K)
    C_v  :: Float64   # Water vapour concentration        (mol·m⁻³)
    s    :: Float64   # Liquid water saturation           (–)
    C_H2 :: Float64   # Hydrogen concentration            (mol·m⁻³)
    C_N2 :: Float64   # Nitrogen concentration            (mol·m⁻³)
    C_O2::Float64   # Additional gas inventory (mol/m^3)
    C_CO2::Float64 # CO₂ gas inventory
end

"""Internal state at one anode microporous-layer location."""
struct AnodeMPLState <: AbstractCellState
    T    :: Float64   # Temperature                       (K)
    C_v  :: Float64   # Water vapour concentration        (mol·m⁻³)
    s    :: Float64   # Liquid water saturation           (–)
    C_H2 :: Float64   # Hydrogen concentration            (mol·m⁻³)
    C_N2 :: Float64   # Nitrogen concentration            (mol·m⁻³)
    C_O2::Float64   # Additional gas inventory (mol/m^3)
    C_CO2::Float64 # CO₂ gas inventory
end

"""Internal state at the anode catalyst-layer location.
Contains ionomer water content and the anode electrode potential."""
struct AnodeCLState <: AbstractCellState
    T      :: Float64   # Temperature                     (K)
    C_v    :: Float64   # Water vapour concentration      (mol·m⁻³)
    s      :: Float64   # Liquid water saturation         (–)
    lambda :: Float64   # Ionomer water content           (–)
    C_H2   :: Float64   # Hydrogen concentration          (mol·m⁻³)
    C_N2   :: Float64   # Nitrogen concentration          (mol·m⁻³)
    C_O2::Float64   # Additional gas inventory (mol/m^3)
    phi_a::Float64  # Anode electrode potential vs SHE (V)
    theta_PtOH::Float64
    theta_Pt_sO::Float64
    theta_Pt_bO::Float64
    C_CO2::Float64 # CO₂ gas inventory
end

# ────────────────────────────────────────────────────────────────────────────────
# Electrolyte
# ────────────────────────────────────────────────────────────────────────────────

"""Internal state at the membrane location.
Only `T` and `lambda` are defined here: gas crossover is represented by fluxes
rather than a membrane gas inventory, and liquid water has no separate phase."""
struct MembraneState <: AbstractCellState
    T      :: Float64   # Temperature                     (K)
    lambda :: Float64   # Ionomer water content           (–)
end

# ────────────────────────────────────────────────────────────────────────────────
# Cathode side
# ────────────────────────────────────────────────────────────────────────────────

"""Internal state at the cathode catalyst-layer location.
Contains ionomer water content, the legacy ORR overpotential and the cathode
electrode potential."""
struct CathodeCLState <: AbstractCellState
    T      :: Float64   # Temperature                     (K)
    C_v    :: Float64   # Water vapour concentration      (mol·m⁻³)
    s      :: Float64   # Liquid water saturation         (–)
    lambda :: Float64   # Ionomer water content           (–)
    C_O2   :: Float64   # Oxygen concentration            (mol·m⁻³)
    C_N2   :: Float64   # Nitrogen concentration          (mol·m⁻³)
    eta_c  :: Float64   # Cathode overpotential           (V)
    C_H2::Float64   # Additional gas inventory (mol/m^3)
    phi_c::Float64  # Cathode electrode potential vs SHE (V)
    theta_PtOH::Float64
    theta_Pt_sO::Float64
    theta_Pt_bO::Float64
    C_CO2::Float64 # CO₂ gas inventory
end

"""Internal state at one cathode microporous-layer location."""
struct CathodeMPLState <: AbstractCellState
    T    :: Float64   # Temperature                       (K)
    C_v  :: Float64   # Water vapour concentration        (mol·m⁻³)
    s    :: Float64   # Liquid water saturation           (–)
    C_O2 :: Float64   # Oxygen concentration              (mol·m⁻³)
    C_N2 :: Float64   # Nitrogen concentration            (mol·m⁻³)
    C_H2::Float64   # Additional gas inventory (mol/m^3)
    C_CO2::Float64 # CO₂ gas inventory
end

"""Internal state at one cathode gas-diffusion-layer location."""
struct CathodeGDLState <: AbstractCellState
    T    :: Float64   # Temperature                       (K)
    C_v  :: Float64   # Water vapour concentration        (mol·m⁻³)
    s    :: Float64   # Liquid water saturation           (–)
    C_O2 :: Float64   # Oxygen concentration              (mol·m⁻³)
    C_N2 :: Float64   # Nitrogen concentration            (mol·m⁻³)
    C_H2::Float64   # Additional gas inventory (mol/m^3)
    C_CO2::Float64 # CO₂ gas inventory
end

"""Internal state at one cathode gas-channel location."""
struct CathodeGCState <: AbstractCellState
    T    :: Float64   # Temperature                       (K)
    C_v  :: Float64   # Water vapour concentration        (mol·m⁻³)
    s    :: Float64   # Liquid water saturation           (–)
    C_O2 :: Float64   # Oxygen concentration              (mol·m⁻³)
    C_N2 :: Float64   # Nitrogen concentration            (mol·m⁻³)
    C_H2::Float64   # Additional gas inventory (mol/m^3)
    C_CO2::Float64 # CO₂ gas inventory
end

# ────────────────────────────────────────────────────────────────────────────────
# Manifold nodes  (mixture state: pressure P, relative humidity Phi)
# ────────────────────────────────────────────────────────────────────────────────

"""Internal state at one supply or exhaust manifold location.
Manifolds are characterised by mixture properties (P, Phi) rather than species concentrations.
"""
struct ManifoldState <: AbstractCellState
    P::Float64     # Pressure                            (Pa)
    Phi::Float64   # Relative humidity                   (–)
end

# ────────────────────────────────────────────────────────────────────────────────
# 1D cell-column state (MEA + AGC/CGC, one column = one GC node)
# ────────────────────────────────────────────────────────────────────────────────

"""Complete 1D internal state for one gas-channel column.

Type parameters
---------------
nb_gdl : Int   Number of nodes in each GDL (anode and cathode share the same count).
nb_mpl : Int   Number of nodes in each MPL.

Fields follow the through-plane order from anode to cathode.
"""
struct CellState1D{nb_gdl, nb_mpl}
    agc  :: AnodeGCState
    agdl :: NTuple{nb_gdl, AnodeGDLState}
    ampl :: NTuple{nb_mpl, AnodeMPLState}
    acl  :: AnodeCLState
    mem  :: MembraneState
    ccl  :: CathodeCLState
    cmpl :: NTuple{nb_mpl, CathodeMPLState}
    cgdl :: NTuple{nb_gdl, CathodeGDLState}
    cgc  :: CathodeGCState
end

# ────────────────────────────────────────────────────────────────────────────────
# Manifold line (one per manifold: anode supply, anode exhaust, cathode supply, cathode exhaust)
# ────────────────────────────────────────────────────────────────────────────────

"""Complete state of one manifold line (supply or exhaust).

Type parameters
---------------
nb_nodes : Int   Number of spatial nodes in this manifold.
                 Currently nb_nodes=1 (single-node manifold), but extensible.
"""
struct ManifoldLine{nb_nodes}
    nodes :: NTuple{nb_nodes, ManifoldState}
end

"""Typed bundle for the four manifold state lines.

This lightweight container keeps manifold entities grouped while preserving
separate manifold lines (asm, aem, csm, cem).
"""
struct _ManifoldStateBundle{ASM, AEM, CSM, CEM}
    asm::ASM
    aem::AEM
    csm::CSM
    cem::CEM
end

"""Typed bundle for the four manifold derivative lines.

Field types are generic to avoid coupling with declaration order between
node and equation model files.
"""
struct _ManifoldDerivativeBundle{ASM, AEM, CSM, CEM}
    asm::ASM
    aem::AEM
    csm::CSM
    cem::CEM
end
# ────────────────────────────────────────────────────────────────────────────────
# P2D fuel-cell state (cell columns = MEA + AGC/CGC)
# ────────────────────────────────────────────────────────────────────────────────

"""Complete P2D (1D+1D) internal state of the fuel-cell cell-column stack.
Each node includes the full through-plane column (AGC + MEA core + CGC).

Type parameters
---------------
nb_gdl : Int   Number of nodes per GDL.
nb_mpl : Int   Number of nodes per MPL.
nb_gc  : Int   Number of gas-channel nodes (spatial discretisation along the channel).
               Typical range: 1 – 10.

Fields:
  - `nodes`: cell-column state at each gas-channel position.
"""
struct FuelCellStateP2D{nb_gdl, nb_mpl, nb_gc}
    nodes :: NTuple{nb_gc, CellState1D{nb_gdl, nb_mpl}}
end



# Existing constructors initialize the additional species to zero.
AnodeGCState(T, C_v, s, C_H2, C_N2) = AnodeGCState(T, C_v, s, C_H2, C_N2, 0.0)

# Existing constructors initialize the additional species to zero.
AnodeGDLState(T, C_v, s, C_H2, C_N2) = AnodeGDLState(T, C_v, s, C_H2, C_N2, 0.0)

# Existing constructors initialize the additional species to zero.
AnodeMPLState(T, C_v, s, C_H2, C_N2) = AnodeMPLState(T, C_v, s, C_H2, C_N2, 0.0)

# Existing constructors initialize the additional species to zero.
AnodeCLState(T, C_v, s, lambda, C_H2, C_N2) = AnodeCLState(T, C_v, s, lambda, C_H2, C_N2, 0.0, 0.0)
AnodeCLState(T, C_v, s, lambda, C_H2, C_N2, C_O2) = AnodeCLState(T, C_v, s, lambda, C_H2, C_N2, C_O2, 0.0)
AnodeCLState(T, C_v, s, lambda, C_H2, C_N2, C_O2, phi_a) =
    AnodeCLState(T, C_v, s, lambda, C_H2, C_N2, C_O2, phi_a, 0.0, 0.0, 0.0)

# Existing constructors initialize the additional species to zero.
CathodeCLState(T, C_v, s, lambda, C_O2, C_N2, eta_c) = CathodeCLState(T, C_v, s, lambda, C_O2, C_N2, eta_c, 0.0, 0.0)
CathodeCLState(T, C_v, s, lambda, C_O2, C_N2, eta_c, C_H2) = CathodeCLState(T, C_v, s, lambda, C_O2, C_N2, eta_c, C_H2, 0.0)
CathodeCLState(T, C_v, s, lambda, C_O2, C_N2, eta_c, C_H2, phi_c) =
    CathodeCLState(T, C_v, s, lambda, C_O2, C_N2, eta_c, C_H2, phi_c, 0.0, 0.0, 0.0)

# Existing constructors initialize the additional species to zero.
CathodeMPLState(T, C_v, s, C_O2, C_N2) = CathodeMPLState(T, C_v, s, C_O2, C_N2, 0.0)

# Existing constructors initialize the additional species to zero.
CathodeGDLState(T, C_v, s, C_O2, C_N2) = CathodeGDLState(T, C_v, s, C_O2, C_N2, 0.0)

# Existing constructors initialize the additional species to zero.
CathodeGCState(T, C_v, s, C_O2, C_N2) = CathodeGCState(T, C_v, s, C_O2, C_N2, 0.0)

# Compatibility: existing constructors initialize CO₂ to zero.
AnodeGCState(x0, x1, x2, x3, x4, x5) = AnodeGCState(x0, x1, x2, x3, x4, x5, 0.0)
AnodeGDLState(x0, x1, x2, x3, x4, x5) = AnodeGDLState(x0, x1, x2, x3, x4, x5, 0.0)
AnodeMPLState(x0, x1, x2, x3, x4, x5) = AnodeMPLState(x0, x1, x2, x3, x4, x5, 0.0)
AnodeCLState(x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10) = AnodeCLState(x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, 0.0)
CathodeCLState(x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11) = CathodeCLState(x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, 0.0)
CathodeMPLState(x0, x1, x2, x3, x4, x5) = CathodeMPLState(x0, x1, x2, x3, x4, x5, 0.0)
CathodeGDLState(x0, x1, x2, x3, x4, x5) = CathodeGDLState(x0, x1, x2, x3, x4, x5, 0.0)
CathodeGCState(x0, x1, x2, x3, x4, x5) = CathodeGCState(x0, x1, x2, x3, x4, x5, 0.0)
