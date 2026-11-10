# -*- coding: utf-8 -*-

"""
    AlphaPEM.Config

This module is the configuration entry point for AlphaPEM,
and provides access to the main configuration parameters.

Modules:
    - current_parameters: Structures for current density parameters
    - fuel_cell_parameters: Structures for physical, operating, numerical, and experimental parameters
    - simulation_config: Structure and validation for simulation configuration parameters
    - calibration_config: Structures for Genetic Algorithm calibration configuration
    - reaction_parameters: Independent local HOR/ORR reproduction parameters

Exports:
    - AbstractCurrentParams, StepParams, PolarizationParams, PolarizationCalibrationParams, EISParams
    - AbstractFuelCellParams, PhysicalParams, OperatingConditions, PolaExperimentalData, NumericalParams
    - SimulationConfig, validate_config
    - GAConfig, CalibrationConfig, CalibrationResult
    - ReactionParams, validate_reaction_parameters
"""
module Config

using ..Utils: y_O2_ext

include("current_parameters.jl")
include("fuel_cell_parameters.jl")
include("numerical_parameters.jl")
include("reaction_parameters.jl")
include("state_scaling.jl")
include("simulation_config.jl")

using .StateScalingModule: CellStateScaling, ManifoldStateScaling, AuxiliaryStateScaling,
                           CurrentDistributionScaling, DAEAlgebraicScaling, StateScaling
using .SimulationConfigModule: SimulationConfig, validate_config

include("calibration_config.jl")

export AbstractCurrentParams, StepParams, PolarizationParams, PolarizationCalibrationParams, EISParams
export validate_gas_feed
export ReactionParams, validate_reaction_parameters
export AbstractFuelCellParams, PhysicalParams, OperatingConditions, PolaExperimentalData, NumericalParams,
       PARAMETER_METADATA, UNDETERMINED_PARAMETER_BOUNDS
export CellStateScaling, ManifoldStateScaling, AuxiliaryStateScaling,
       CurrentDistributionScaling, DAEAlgebraicScaling, StateScaling
export SimulationConfig, validate_config
export GAConfig, CalibrationConfig, CalibrationResult

end  # module Config

