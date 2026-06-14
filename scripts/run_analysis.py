from __future__ import annotations

import argparse
import csv
import os
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

os.environ.setdefault("MPLBACKEND", "Agg")
os.environ.setdefault("MPLCONFIGDIR", str(ROOT / "results" / ".mplconfig"))

import matplotlib.pyplot as plt

from rocket_combustion.model import BAR, run_of_sweep, run_pressure_sweep


def write_csv(path: Path, rows: list[dict[str, float]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def column(rows: list[dict[str, float]], name: str) -> np.ndarray:
    return np.array([row[name] for row in rows], dtype=float)


def save_temperature_plot(rows: list[dict[str, float]], figures_dir: Path) -> None:
    of = column(rows, "of_ratio")
    temperature = column(rows, "chamber_temperature")

    fig, ax = plt.subplots(figsize=(7.0, 4.4))
    ax.plot(of, temperature, marker="o", markersize=3.0, linewidth=1.8)
    ax.axvline(4.0, color="0.35", linestyle="--", linewidth=1.0, label="stoichiometric O/F")
    ax.set_xlabel("Oxidizer-to-fuel mass ratio, O/F [-]")
    ax.set_ylabel("Adiabatic chamber temperature [K]")
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(figures_dir / "chamber_temperature_vs_of.png", dpi=220)
    plt.close(fig)


def save_performance_plot(rows: list[dict[str, float]], figures_dir: Path) -> None:
    of = column(rows, "of_ratio")
    isp_sl = column(rows, "specific_impulse")
    isp_vac = column(rows, "specific_impulse_vacuum")
    cstar = column(rows, "characteristic_velocity")

    fig, ax1 = plt.subplots(figsize=(7.0, 4.4))
    ax1.plot(of, isp_sl, color="#0b6e4f", marker="o", markersize=3.0, label="sea-level Isp")
    ax1.plot(of, isp_vac, color="#2563eb", marker="^", markersize=3.0, label="vacuum Isp")
    ax1.set_xlabel("Oxidizer-to-fuel mass ratio, O/F [-]")
    ax1.set_ylabel("Ideal specific impulse [s]")
    ax1.grid(True, alpha=0.3)

    ax2 = ax1.twinx()
    ax2.plot(of, cstar, color="#9a3412", marker="s", markersize=2.8, label="c*")
    ax2.set_ylabel("Characteristic velocity c* [m/s]", color="#9a3412")
    ax2.tick_params(axis="y", labelcolor="#9a3412")

    lines = ax1.get_lines() + ax2.get_lines()
    ax1.legend(lines, [line.get_label() for line in lines], loc="best")
    fig.tight_layout()
    fig.savefig(figures_dir / "performance_vs_of.png", dpi=220)
    plt.close(fig)


def save_gamma_mw_plot(rows: list[dict[str, float]], figures_dir: Path) -> None:
    of = column(rows, "of_ratio")
    gamma = column(rows, "gamma")
    mw = column(rows, "molecular_weight")

    fig, ax1 = plt.subplots(figsize=(7.0, 4.4))
    ax1.plot(of, gamma, color="#1d4ed8", marker="o", markersize=3.0, label="gamma")
    ax1.set_xlabel("Oxidizer-to-fuel mass ratio, O/F [-]")
    ax1.set_ylabel("Heat capacity ratio gamma [-]", color="#1d4ed8")
    ax1.tick_params(axis="y", labelcolor="#1d4ed8")
    ax1.grid(True, alpha=0.3)

    ax2 = ax1.twinx()
    ax2.plot(of, mw, color="#7c3aed", marker="s", markersize=2.8, label="molecular weight")
    ax2.set_ylabel("Mean molecular weight [kg/kmol]", color="#7c3aed")
    ax2.tick_params(axis="y", labelcolor="#7c3aed")

    lines = ax1.get_lines() + ax2.get_lines()
    ax1.legend(lines, [line.get_label() for line in lines], loc="best")
    fig.tight_layout()
    fig.savefig(figures_dir / "gamma_and_molecular_weight_vs_of.png", dpi=220)
    plt.close(fig)


def save_species_plot(rows: list[dict[str, float]], figures_dir: Path) -> None:
    of = column(rows, "of_ratio")
    species = {
        "CO2": "mole_fraction_co2",
        "H2O": "mole_fraction_h2o",
        "CO": "mole_fraction_co",
        "H2": "mole_fraction_h2",
        "O2": "mole_fraction_o2",
    }

    fig, ax = plt.subplots(figsize=(7.0, 4.4))
    for label, field in species.items():
        ax.plot(of, column(rows, field), marker="o", markersize=2.5, linewidth=1.5, label=label)
    ax.set_xlabel("Oxidizer-to-fuel mass ratio, O/F [-]")
    ax.set_ylabel("Equilibrium mole fraction [-]")
    ax.grid(True, alpha=0.3)
    ax.legend(ncol=3)
    fig.tight_layout()
    fig.savefig(figures_dir / "major_species_vs_of.png", dpi=220)
    plt.close(fig)


def save_pressure_plot(rows: list[dict[str, float]], figures_dir: Path) -> None:
    pressure = column(rows, "chamber_pressure_bar")
    temperature = column(rows, "chamber_temperature")
    isp_sl = column(rows, "specific_impulse")
    isp_vac = column(rows, "specific_impulse_vacuum")

    fig, ax1 = plt.subplots(figsize=(7.0, 4.4))
    ax1.plot(pressure, temperature, color="#b45309", marker="o", label="temperature")
    ax1.set_xlabel("Chamber pressure [bar]")
    ax1.set_ylabel("Adiabatic chamber temperature [K]", color="#b45309")
    ax1.tick_params(axis="y", labelcolor="#b45309")
    ax1.grid(True, alpha=0.3)

    ax2 = ax1.twinx()
    ax2.plot(pressure, isp_sl, color="#0f766e", marker="s", label="sea-level Isp")
    ax2.plot(pressure, isp_vac, color="#2563eb", marker="^", label="vacuum Isp")
    ax2.set_ylabel("Ideal specific impulse [s]", color="#0f766e")
    ax2.tick_params(axis="y", labelcolor="#0f766e")

    lines = ax1.get_lines() + ax2.get_lines()
    ax1.legend(lines, [line.get_label() for line in lines], loc="best")
    fig.tight_layout()
    fig.savefig(figures_dir / "pressure_sweep.png", dpi=220)
    plt.close(fig)


def save_transport_plot(rows: list[dict[str, float]], figures_dir: Path) -> None:
    of = column(rows, "of_ratio")
    viscosity = column(rows, "viscosity") * 1e6
    prandtl = column(rows, "prandtl")

    fig, ax1 = plt.subplots(figsize=(7.0, 4.4))
    ax1.plot(of, viscosity, color="#be123c", marker="o", markersize=3.0, label="dynamic viscosity")
    ax1.set_xlabel("Oxidizer-to-fuel mass ratio, O/F [-]")
    ax1.set_ylabel("Dynamic viscosity [uPa s]", color="#be123c")
    ax1.tick_params(axis="y", labelcolor="#be123c")
    ax1.grid(True, alpha=0.3)

    ax2 = ax1.twinx()
    ax2.plot(of, prandtl, color="#4338ca", marker="s", markersize=2.8, label="Prandtl number")
    ax2.set_ylabel("Prandtl number [-]", color="#4338ca")
    ax2.tick_params(axis="y", labelcolor="#4338ca")

    lines = ax1.get_lines() + ax2.get_lines()
    ax1.legend(lines, [line.get_label() for line in lines], loc="best")
    fig.tight_layout()
    fig.savefig(figures_dir / "transport_vs_of.png", dpi=220)
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Run LOX/LCH4 equilibrium rocket-combustion analysis."
    )
    parser.add_argument("--output-dir", type=Path, default=ROOT / "results")
    parser.add_argument("--chamber-pressure-bar", type=float, default=20.0)
    parser.add_argument("--area-ratio", type=float, default=15.0)
    parser.add_argument("--ambient-pressure-bar", type=float, default=1.0)
    args = parser.parse_args()

    output_dir = args.output_dir
    figures_dir = output_dir / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True)

    of_ratios = np.linspace(2.0, 4.4, 49)
    of_rows = run_of_sweep(
        of_ratios,
        chamber_pressure=args.chamber_pressure_bar * BAR,
        nozzle_area_ratio=args.area_ratio,
        ambient_pressure=args.ambient_pressure_bar * BAR,
    )
    write_csv(output_dir / "of_sweep.csv", of_rows)

    best_row = max(of_rows, key=lambda row: row["specific_impulse_vacuum"])
    best_of = best_row["of_ratio"]
    pressure_values = np.array([20.0, 40.0, 60.0, 80.0, 100.0]) * BAR
    pressure_rows = run_pressure_sweep(
        pressure_values,
        of_ratio=best_of,
        nozzle_area_ratio=args.area_ratio,
        ambient_pressure=args.ambient_pressure_bar * BAR,
    )
    write_csv(output_dir / "pressure_sweep.csv", pressure_rows)

    save_temperature_plot(of_rows, figures_dir)
    save_performance_plot(of_rows, figures_dir)
    save_gamma_mw_plot(of_rows, figures_dir)
    save_species_plot(of_rows, figures_dir)
    save_pressure_plot(pressure_rows, figures_dir)
    save_transport_plot(of_rows, figures_dir)

    print("Analysis complete.")
    print(f"Best ideal vacuum Isp: {best_row['specific_impulse_vacuum']:.2f} s at O/F={best_of:.2f}")
    print(f"Sea-level Isp at that point: {best_row['specific_impulse']:.2f} s")
    print(f"Maximum chamber temperature: {max(column(of_rows, 'chamber_temperature')):.1f} K")
    print(f"Results written to: {output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
