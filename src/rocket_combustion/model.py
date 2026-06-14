from __future__ import annotations

from dataclasses import asdict, dataclass
from math import isfinite, sqrt
from typing import Iterable

import cantera as ct


G0 = 9.80665
BAR = 100_000.0


@dataclass(frozen=True)
class RocketCase:
    """Input parameters for one equilibrium rocket-combustion case."""

    of_ratio: float
    chamber_pressure: float = 20.0 * BAR
    inlet_temperature: float = 300.0
    fuel: str = "CH4"
    oxidizer: str = "O2"
    mechanism: str = "gri30_highT.yaml"
    nozzle_area_ratio: float = 15.0
    ambient_pressure: float = 1.0 * BAR


@dataclass(frozen=True)
class ChamberState:
    of_ratio: float
    chamber_pressure_bar: float
    chamber_temperature: float
    gamma: float
    molecular_weight: float
    gas_constant: float
    density: float
    cp_mass: float
    cv_mass: float
    sound_speed: float
    viscosity: float
    thermal_conductivity: float
    prandtl: float
    mole_fraction_co2: float
    mole_fraction_h2o: float
    mole_fraction_co: float
    mole_fraction_h2: float
    mole_fraction_o2: float


@dataclass(frozen=True)
class NozzlePerformance:
    exit_mach: float
    exit_pressure_bar: float
    exit_temperature: float
    exit_velocity: float
    thrust_coefficient: float
    thrust_coefficient_vacuum: float
    characteristic_velocity: float
    specific_impulse: float
    specific_impulse_vacuum: float
    ambient_pressure_bar: float
    nozzle_area_ratio: float


def mole_fractions_from_of(
    gas: ct.Solution, fuel: str, oxidizer: str, of_ratio: float
) -> dict[str, float]:
    """Return reactant mole fractions corresponding to a mass O/F ratio."""

    if of_ratio <= 0.0:
        raise ValueError("O/F ratio must be positive.")

    fuel_index = gas.species_index(fuel)
    oxidizer_index = gas.species_index(oxidizer)
    fuel_mw = gas.molecular_weights[fuel_index]
    oxidizer_mw = gas.molecular_weights[oxidizer_index]

    fuel_mass = 1.0
    oxidizer_mass = of_ratio * fuel_mass
    fuel_kmol = fuel_mass / fuel_mw
    oxidizer_kmol = oxidizer_mass / oxidizer_mw
    total_kmol = fuel_kmol + oxidizer_kmol

    return {
        fuel: fuel_kmol / total_kmol,
        oxidizer: oxidizer_kmol / total_kmol,
    }


def _species_mole_fraction(gas: ct.Solution, species: str) -> float:
    try:
        return float(gas.X[gas.species_index(species)])
    except ValueError:
        return 0.0


def equilibrate_chamber(case: RocketCase) -> tuple[ct.Solution, ChamberState]:
    """Solve an adiabatic, constant-pressure combustion state."""

    gas = ct.Solution(case.mechanism)
    composition = mole_fractions_from_of(
        gas, case.fuel, case.oxidizer, case.of_ratio
    )
    gas.TPX = case.inlet_temperature, case.chamber_pressure, composition
    gas.equilibrate("HP")

    gamma = gas.cp_mass / gas.cv_mass
    gas_constant = ct.gas_constant / gas.mean_molecular_weight
    sound_speed = sqrt(gamma * gas_constant * gas.T)
    thermal_conductivity = gas.thermal_conductivity
    viscosity = gas.viscosity

    state = ChamberState(
        of_ratio=case.of_ratio,
        chamber_pressure_bar=case.chamber_pressure / BAR,
        chamber_temperature=gas.T,
        gamma=gamma,
        molecular_weight=gas.mean_molecular_weight,
        gas_constant=gas_constant,
        density=gas.density,
        cp_mass=gas.cp_mass,
        cv_mass=gas.cv_mass,
        sound_speed=sound_speed,
        viscosity=viscosity,
        thermal_conductivity=thermal_conductivity,
        prandtl=gas.cp_mass * viscosity / thermal_conductivity,
        mole_fraction_co2=_species_mole_fraction(gas, "CO2"),
        mole_fraction_h2o=_species_mole_fraction(gas, "H2O"),
        mole_fraction_co=_species_mole_fraction(gas, "CO"),
        mole_fraction_h2=_species_mole_fraction(gas, "H2"),
        mole_fraction_o2=_species_mole_fraction(gas, "O2"),
    )
    return gas, state


def area_mach_relation(mach: float, gamma: float) -> float:
    """Isentropic area ratio A/A* for a calorically perfect gas."""

    if mach <= 0.0:
        raise ValueError("Mach number must be positive.")
    exponent = (gamma + 1.0) / (2.0 * (gamma - 1.0))
    bracket = (2.0 / (gamma + 1.0)) * (1.0 + 0.5 * (gamma - 1.0) * mach**2)
    return (1.0 / mach) * bracket**exponent


def supersonic_mach_from_area_ratio(area_ratio: float, gamma: float) -> float:
    """Find the supersonic Mach number for a specified nozzle area ratio."""

    if area_ratio < 1.0:
        raise ValueError("Area ratio must be greater than or equal to one.")

    low = 1.0 + 1e-9
    high = 2.0
    while area_mach_relation(high, gamma) < area_ratio:
        high *= 1.5
        if high > 100.0:
            raise RuntimeError("Could not bracket supersonic Mach number.")

    for _ in range(100):
        mid = 0.5 * (low + high)
        if area_mach_relation(mid, gamma) < area_ratio:
            low = mid
        else:
            high = mid
    return 0.5 * (low + high)


def ideal_characteristic_velocity(temperature: float, gamma: float, gas_constant: float) -> float:
    """Ideal frozen-gamma characteristic velocity."""

    return sqrt(gas_constant * temperature / gamma) * (
        (gamma + 1.0) / 2.0
    ) ** ((gamma + 1.0) / (2.0 * (gamma - 1.0)))


def ideal_nozzle_performance(case: RocketCase, chamber: ChamberState) -> NozzlePerformance:
    """Compute ideal frozen-gamma nozzle performance for a fixed area ratio."""

    gamma = chamber.gamma
    mach_exit = supersonic_mach_from_area_ratio(case.nozzle_area_ratio, gamma)
    pressure_ratio = (
        1.0 + 0.5 * (gamma - 1.0) * mach_exit**2
    ) ** (-gamma / (gamma - 1.0))
    temperature_ratio = (
        1.0 + 0.5 * (gamma - 1.0) * mach_exit**2
    ) ** -1.0

    exit_pressure = case.chamber_pressure * pressure_ratio
    exit_temperature = chamber.chamber_temperature * temperature_ratio
    exit_velocity = mach_exit * sqrt(gamma * chamber.gas_constant * exit_temperature)
    pressure_thrust = (
        (exit_pressure - case.ambient_pressure)
        * case.nozzle_area_ratio
        / case.chamber_pressure
    )
    pressure_thrust_vacuum = exit_pressure * case.nozzle_area_ratio / case.chamber_pressure
    momentum_thrust = sqrt(
        (2.0 * gamma**2 / (gamma - 1.0))
        * (2.0 / (gamma + 1.0)) ** ((gamma + 1.0) / (gamma - 1.0))
        * (1.0 - pressure_ratio ** ((gamma - 1.0) / gamma))
    )
    thrust_coefficient = momentum_thrust + pressure_thrust
    thrust_coefficient_vacuum = momentum_thrust + pressure_thrust_vacuum
    characteristic_velocity = ideal_characteristic_velocity(
        chamber.chamber_temperature, gamma, chamber.gas_constant
    )

    return NozzlePerformance(
        exit_mach=mach_exit,
        exit_pressure_bar=exit_pressure / BAR,
        exit_temperature=exit_temperature,
        exit_velocity=exit_velocity,
        thrust_coefficient=thrust_coefficient,
        thrust_coefficient_vacuum=thrust_coefficient_vacuum,
        characteristic_velocity=characteristic_velocity,
        specific_impulse=characteristic_velocity * thrust_coefficient / G0,
        specific_impulse_vacuum=characteristic_velocity
        * thrust_coefficient_vacuum
        / G0,
        ambient_pressure_bar=case.ambient_pressure / BAR,
        nozzle_area_ratio=case.nozzle_area_ratio,
    )


def run_case(case: RocketCase) -> dict[str, float]:
    """Run one chamber/nozzle case and return a flat results dictionary."""

    _, chamber = equilibrate_chamber(case)
    nozzle = ideal_nozzle_performance(case, chamber)
    row = asdict(chamber)
    row.update(asdict(nozzle))
    return row


def run_of_sweep(
    of_ratios: Iterable[float],
    chamber_pressure: float = 20.0 * BAR,
    nozzle_area_ratio: float = 15.0,
    ambient_pressure: float = 1.0 * BAR,
) -> list[dict[str, float]]:
    """Run a sweep over oxidizer-to-fuel mass ratio."""

    rows = []
    for of_ratio in of_ratios:
        case = RocketCase(
            of_ratio=float(of_ratio),
            chamber_pressure=chamber_pressure,
            nozzle_area_ratio=nozzle_area_ratio,
            ambient_pressure=ambient_pressure,
        )
        row = run_case(case)
        if not all(isfinite(value) for value in row.values()):
            raise FloatingPointError(f"Non-finite result for O/F={of_ratio}.")
        rows.append(row)
    return rows


def run_pressure_sweep(
    pressures: Iterable[float],
    of_ratio: float,
    nozzle_area_ratio: float = 15.0,
    ambient_pressure: float = 1.0 * BAR,
) -> list[dict[str, float]]:
    """Run a sweep over chamber pressure for a selected O/F ratio."""

    rows = []
    for pressure in pressures:
        case = RocketCase(
            of_ratio=of_ratio,
            chamber_pressure=float(pressure),
            nozzle_area_ratio=nozzle_area_ratio,
            ambient_pressure=ambient_pressure,
        )
        row = run_case(case)
        if not all(isfinite(value) for value in row.values()):
            raise FloatingPointError(
                f"Non-finite result for pressure={pressure} Pa."
            )
        rows.append(row)
    return rows
