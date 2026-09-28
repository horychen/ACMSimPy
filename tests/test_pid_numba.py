"""Check nopython execution using the tutorial's actual regulator jitclass."""

import ast
import unittest

from test_pid_regression import ROOT, controllers

try:
    from numba import float64, njit
    from numba.experimental import jitclass
except ImportError:
    njit = None


@unittest.skipIf(njit is None, "Numba is not installed")
class PidNumbaRegression(unittest.TestCase):
    def test_compiled_controllers_match_python(self):
        path = ROOT / "simulation" / "tutorials_ep5_pid_regulation.py"
        tree = ast.parse(path.read_text(encoding="utf-8"))
        node = next(node for node in tree.body
                    if isinstance(node, ast.ClassDef) and node.name == "The_PID_Regulator")
        namespace = {"jitclass": jitclass, "float64": float64}
        exec(compile(ast.Module(body=[node], type_ignores=[]), str(path), "exec"), namespace)
        regulator = namespace["The_PID_Regulator"]

        for name, pid in controllers("tustin_pid"):
            with self.subTest(file=name):
                compiled = njit(pid)
                reference = regulator(1.0, 2.0, 0.5, 0.1, 10.0, 2.0, 0.01)
                actual = regulator(1.0, 2.0, 0.5, 0.1, 10.0, 2.0, 0.01)
                for setpoint, measurement in ((30.0, 0.0), (0.0, 1.0), (0.0, 1.0),
                                              (-30.0, -1.0), (0.0, 0.0)):
                    reference.setpoint = actual.setpoint = setpoint
                    reference.measurement = actual.measurement = measurement
                    self.assertAlmostEqual(compiled(actual), pid(reference))
                    self.assertAlmostEqual(actual.integrator, reference.integrator)
                    self.assertAlmostEqual(actual.differentiator, reference.differentiator)
                self.assertTrue(compiled.nopython_signatures)


if __name__ == "__main__":
    unittest.main()
