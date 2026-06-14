# Equilibrium Performance Analysis of a LOX/LCH4 Rocket Engine Using Cantera

Project for the course **Computer Methods in Combustion**.

This repository contains a Python/Cantera model of an idealized methane-oxygen rocket combustion chamber. The project investigates how oxidizer-to-fuel mass ratio and chamber pressure influence combustion temperature, equilibrium product composition, gas properties and ideal nozzle performance.

The calculations are intended as a preliminary engineering analysis, not as a complete rocket engine design tool.

## Project Scope

The model calculates:

- adiabatic equilibrium chamber temperature,
- equilibrium mole fractions of major combustion products,
- heat capacity ratio, mean molecular weight and gas constant,
- density, sound speed and selected transport properties,
- ideal characteristic velocity `c*`,
- ideal thrust coefficient,
- sea-level and vacuum specific impulse.

The baseline propellant combination is:

- fuel: methane, `CH4`,
- oxidizer: oxygen, `O2`,
- mechanism: Cantera `gri30_highT.yaml`,
- baseline chamber pressure: 20 bar,
- baseline nozzle area ratio: 15.



## Installation

Create a virtual environment and install the required packages:

```powershell
python -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r requirements.txt
```

## Running the Analysis

Run the main analysis script:

```powershell
.\.venv\Scripts\python.exe scripts\run_analysis.py
```

The script generates:

- `results/of_sweep.csv`,
- `results/pressure_sweep.csv`,
- plots in `results/figures/`.

## Running Tests

```powershell
.\.venv\Scripts\python.exe -m unittest discover -s tests
```

The tests check:

- O/F ratio conversion,
- isentropic area-Mach relation inversion,
- basic physical consistency of a representative chamber and nozzle case.

## Report

The final report is included as:

```text
MKWS_2026_Michna_R_raport.pdf
```

## Main Findings

For the baseline case used in the generated results:

- maximum chamber temperature is approximately 3456 K at O/F = 3.70,
- maximum ideal vacuum specific impulse is approximately 337 s at O/F = 2.55,
- maximum temperature and maximum performance occur at different mixture ratios, because molecular weight and expansion properties also affect rocket performance.

## Limitations

The model uses several simplifying assumptions:

- reactants are treated as ideal gases at 300 K,
- the chamber is adiabatic and perfectly mixed,
- only gas-phase chemistry from `gri30_highT.yaml` is included,
- nozzle expansion uses a frozen-gamma ideal-gas approximation,
- no injector losses, heat losses, boundary-layer losses or nozzle efficiencies are included.

Because of these assumptions, the results should be interpreted as ideal trends for coursework and early-stage engineering reasoning.
