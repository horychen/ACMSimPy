"""Exercise tutorial controllers without running their plotting/simulation setup.

Run with: python -m unittest discover -s tests -v
"""

import ast
from pathlib import Path
from types import SimpleNamespace
import unittest


ROOT = Path(__file__).resolve().parents[1]
DYNAMIC = {
    "pid_cjh_taste.py",
    "tutorials_ep9_flux_estimator.py",
    "tutorials_ep10_ParameterIdentification.py",
    "tutorials_ep11_SynIFO.py",
}


def controllers(name):
    for path in sorted((ROOT / "simulation").glob("*.py")):
        source = path.read_text(encoding="utf-8")
        if f"def {name}(" not in source:
            continue
        tree = ast.parse(source, filename=str(path))
        for node in tree.body:
            if isinstance(node, ast.FunctionDef) and node.name == name:
                node.decorator_list = []
                namespace = {}
                exec(compile(ast.Module(body=[node], type_ignores=[]), str(path), "exec"), namespace)
                yield path.name, namespace[name]


def regulator(**changes):
    state = dict(
        Kp=0.0, Ki=0.0, Kd=0.0, tau=1.0, T=1.0,
        OutLimit=100.0, IntLimit=100.0, integrator=0.0,
        prevError=0.0, differentiator=0.0, prevMeasurement=0.0,
        Out=0.0, setpoint=0.0, measurement=0.0,
        OutPrev=0.0, Err=0.0, ErrPrev=0.0, Ref=0.0, Fbk=0.0,
    )
    state.update(changes)
    return SimpleNamespace(**state)


class PidRegression(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.pids = dict(controllers("tustin_pid"))
        cls.pis = dict(controllers("incremental_pi"))

    def test_controller_coverage(self):
        self.assertEqual(len(self.pids), 8)
        self.assertTrue(DYNAMIC <= self.pids.keys())
        self.assertTrue(self.pis)

    def test_derivative_step_matches_bilinear_pole(self):
        for name, pid in self.pids.items():
            for tau in (0.0, 0.25, 0.5, 1.0, 10.0):
                with self.subTest(file=name, tau=tau):
                    reg = regulator(Kd=1.0, tau=tau, measurement=1.0)
                    pole = (2.0 * tau - reg.T) / (2.0 * tau + reg.T)
                    amplitude = -2.0 / (2.0 * tau + reg.T)
                    for k in range(6):
                        pid(reg)
                        self.assertAlmostEqual(reg.differentiator, amplitude * pole**k)

    def test_setpoint_step_does_not_kick_derivative(self):
        for name, pid in self.pids.items():
            with self.subTest(file=name):
                reg = regulator(Kd=1.0, setpoint=10.0)
                pid(reg)
                self.assertEqual(reg.differentiator, 0.0)

    def test_trapezoidal_integral_and_history(self):
        for name, pid in self.pids.items():
            with self.subTest(file=name):
                reg = regulator(Ki=2.0, T=0.1, setpoint=3.0, measurement=1.0, prevError=1.0)
                pid(reg)
                self.assertAlmostEqual(reg.integrator, 0.3)
                self.assertEqual(reg.prevError, 2.0)
                self.assertEqual(reg.prevMeasurement, 1.0)
                pid(reg)
                self.assertAlmostEqual(reg.integrator, 0.7)

    def test_saturation_does_not_create_opposing_integral(self):
        for name in DYNAMIC:
            for sign in (-1.0, 1.0):
                with self.subTest(file=name, sign=sign):
                    reg = regulator(Kp=1.0, Ki=1.0, T=0.1,
                                    setpoint=12.0 * sign, OutLimit=10.0)
                    for _ in range(20):
                        self.pids[name](reg)
                        self.assertEqual(reg.integrator, 0.0)
                        self.assertEqual(reg.Out, 10.0 * sign)
                    # Reducing P releases saturation without an artificial offset.
                    reg.setpoint = 3.0 * sign
                    self.pids[name](reg)
                    self.assertAlmostEqual(reg.integrator, 0.75 * sign)
                    self.assertAlmostEqual(reg.Out, 3.75 * sign)

    def test_saturation_preserves_existing_integral_of_either_sign(self):
        for name in DYNAMIC:
            for sign in (-1.0, 1.0):
                for initial in (-1.0, 1.0):
                    with self.subTest(file=name, sign=sign, initial=initial):
                        reg = regulator(Kp=1.0, Ki=1.0, T=0.1,
                                        setpoint=12.0 * sign, integrator=initial, OutLimit=10.0)
                        self.pids[name](reg)
                        self.assertEqual(reg.integrator, initial)

    def test_integration_can_unwind_while_output_remains_saturated(self):
        for name in DYNAMIC:
            for sign in (-1.0, 1.0):
                with self.subTest(file=name, sign=sign):
                    reg = regulator(Kp=1.0, Ki=1.0, T=0.1, setpoint=-sign,
                                    prevError=-sign, integrator=20.0 * sign, OutLimit=10.0)
                    self.pids[name](reg)
                    self.assertAlmostEqual(reg.integrator, 19.9 * sign)
                    self.assertEqual(reg.Out, 10.0 * sign)

    def test_saturation_check_includes_current_derivative(self):
        for name in DYNAMIC:
            for sign in (-1.0, 1.0):
                with self.subTest(file=name, sign=sign):
                    reg = regulator(Kp=1.0, Ki=1.0, Kd=15.0, T=1.0,
                                    measurement=-sign, setpoint=0.0, OutLimit=10.0)
                    self.pids[name](reg)
                    self.assertEqual(reg.differentiator, 10.0 * sign)
                    self.assertEqual(reg.integrator, 0.0)
                    self.assertEqual(reg.Out, 10.0 * sign)

    def test_tustin_increment_direction_controls_freezing(self):
        for name in DYNAMIC:
            with self.subTest(file=name):
                # Error has reversed, but the trapezoidal increment is still positive.
                reg = regulator(Ki=1.0, T=0.1, setpoint=-1.0, prevError=3.0,
                                integrator=12.0, OutLimit=10.0)
                self.pids[name](reg)
                self.assertEqual(reg.integrator, 12.0)
                self.pids[name](reg)
                self.assertAlmostEqual(reg.integrator, 11.9)

    def test_zero_ki_has_no_integral_memory_or_saturation_offset(self):
        for name, pid in self.pids.items():
            for sign in (-1.0, 1.0):
                with self.subTest(file=name, sign=sign):
                    reg = regulator(Kp=1.0, Ki=1.0, T=0.1, setpoint=sign)
                    pid(reg)
                    self.assertNotEqual(reg.integrator, 0.0)
                    reg.Ki = 0.0
                    for p in (12.0 * sign, 0.0, -12.0 * sign, 0.0):
                        reg.setpoint = p
                        reg.OutLimit = 10.0
                        pid(reg)
                        self.assertEqual(reg.integrator, 0.0)
                        self.assertEqual(reg.Out, max(-10.0, min(10.0, p)))
                    reg.Ki = 1.0
                    reg.setpoint = sign
                    pid(reg)
                    self.assertAlmostEqual(reg.integrator, 0.05 * sign)

    def test_zero_ki_preserves_pd_response(self):
        for name, pid in self.pids.items():
            with self.subTest(file=name):
                reg = regulator(Kp=1.0, Kd=3.0, setpoint=4.0,
                                measurement=1.0, integrator=7.0)
                pid(reg)
                self.assertEqual(reg.integrator, 0.0)
                self.assertAlmostEqual(reg.differentiator, -2.0)
                self.assertAlmostEqual(reg.Out, 1.0)

    def test_first_order_loop_recovers_from_both_saturation_directions(self):
        for name in DYNAMIC:
            for sign in (-1.0, 1.0):
                with self.subTest(file=name, sign=sign):
                    reg = regulator(Kp=2.0, Ki=4.0, T=0.01, OutLimit=1.0)
                    measurement = 0.0
                    # Euler simulation of 0.5 * dy/dt = -y + u.
                    for reference, steps in ((2.0 * sign, 300), (0.25 * sign, 1000)):
                        reg.setpoint = reference
                        for _ in range(steps):
                            reg.measurement = measurement
                            self.pids[name](reg)
                            self.assertLessEqual(abs(reg.Out), 1.0)
                            measurement += reg.T / 0.5 * (reg.Out - measurement)
                    self.assertAlmostEqual(measurement, 0.25 * sign, delta=0.001)
                    self.assertAlmostEqual(reg.integrator, 0.25 * sign, delta=0.001)

    def test_fixed_integral_limits_remain_fixed(self):
        for name in self.pids.keys() - DYNAMIC:
            for initial in (-100.0, 100.0):
                with self.subTest(file=name, initial=initial):
                    reg = regulator(Kp=1.0, Ki=1.0, setpoint=3.0, integrator=initial,
                                    OutLimit=10.0, IntLimit=2.0)
                    self.pids[name](reg)
                    self.assertEqual(reg.integrator, -2.0 if initial < 0.0 else 2.0)
                    self.assertEqual(reg.IntLimit, 2.0)

    def test_incremental_pi_saves_saturated_output_and_recovers(self):
        for name, pi in self.pis.items():
            for sign in (-1.0, 1.0):
                with self.subTest(file=name, sign=sign):
                    reg = regulator(Ki=1.0, setpoint=20.0 * sign, Ref=20.0 * sign, OutLimit=10.0)
                    pi(reg)
                    self.assertEqual(reg.OutPrev, 10.0 * sign)
                    reg.setpoint = -sign
                    reg.Ref = -sign
                    pi(reg)
                    self.assertEqual(reg.Out, 9.0 * sign)
                    self.assertEqual(reg.OutPrev, reg.Out)


if __name__ == "__main__":
    unittest.main()
