import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

import cantera as ct

from rocket_combustion.model import (
    RocketCase,
    area_mach_relation,
    equilibrate_chamber,
    ideal_nozzle_performance,
    mole_fractions_from_of,
    supersonic_mach_from_area_ratio,
)


class RocketCombustionModelTest(unittest.TestCase):
    def test_of_ratio_conversion_is_stoichiometric_for_methane(self):
        gas = ct.Solution("gri30_highT.yaml")
        stoichiometric_of = (
            2.0
            * gas.molecular_weights[gas.species_index("O2")]
            / gas.molecular_weights[gas.species_index("CH4")]
        )
        composition = mole_fractions_from_of(gas, "CH4", "O2", stoichiometric_of)
        self.assertAlmostEqual(composition["O2"] / composition["CH4"], 2.0, places=3)

    def test_area_mach_inverse_on_supersonic_branch(self):
        gamma = 1.22
        mach = supersonic_mach_from_area_ratio(15.0, gamma)
        self.assertGreater(mach, 1.0)
        self.assertAlmostEqual(area_mach_relation(mach, gamma), 15.0, places=6)

    def test_chamber_and_nozzle_are_physical(self):
        case = RocketCase(of_ratio=3.4)
        _, chamber = equilibrate_chamber(case)
        nozzle = ideal_nozzle_performance(case, chamber)

        self.assertGreater(chamber.chamber_temperature, 2500.0)
        self.assertTrue(1.0 < chamber.gamma < 1.4)
        self.assertGreater(nozzle.characteristic_velocity, 1000.0)
        self.assertGreater(nozzle.specific_impulse, 100.0)


if __name__ == "__main__":
    unittest.main()
