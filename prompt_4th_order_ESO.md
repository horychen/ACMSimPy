# Task: Implement 4th-Order ESO with Ramp Load Torque Estimation and Feedforward Compensation

## Goal

Modify the existing PMSM FOC simulation codebase to:
1. Enable the **4th-order position observer** (already partially scaffolded but unused)
2. Wire the observer's disturbance estimate into the speed controller as **feedforward compensation**
3. Apply a **ramp (linearly increasing) load torque** to demonstrate that the 3rd-order observer fails to track ramp disturbance, while the 4th-order observer succeeds
4. Plot comparison figures: 3rd-order vs 4th-order observer under ramp load

## Background: Why 4th-Order?

The mechanical equation of a PMSM is:

```
J * d(omega)/dt = T_em - T_load - B*omega
```

A speed/position observer estimates `omega` and the "total disturbance" `d_to` (which lumps load torque, friction, parameter errors). The observer uses measured rotor position `theta_m` as input.

- **3rd-order observer**: models disturbance as constant (`d(d_to)/dt = 0`). States: `[theta, omega, d_to]`. Can track step load but NOT ramp load.
- **4th-order observer**: models disturbance as ramp (`d(d_to)/dt = p_to`, `d(p_to)/dt = 0`). States: `[theta, omega, d_to, p_to]`. Can track ramp load because it estimates the rate of change of disturbance.

The observer dynamics for 4th-order are:

```
d(x0)/dt = ell1 * e_theta + x1            # x0 = estimated theta
d(x1)/dt = ell2 * e_theta + (Tem + x2) * npp / Js   # x1 = estimated omega (elec.)
d(x2)/dt = ell3 * e_theta + x3            # x2 = estimated disturbance (d_to)
d(x3)/dt = ell4 * e_theta                 # x3 = estimated rate of disturbance (p_to)
```

where `e_theta = theta_measured - x0` is the position estimation error.

Observer gains are placed at `(s + omega_ob)^4`:

```
ell1 = 4 * omega_ob
ell2 = 6 * omega_ob^2
ell3 = 4 * omega_ob^3 * Js / npp
ell4 = omega_ob^4
```

For 3rd-order (current default), gains are at `(s + omega_ob)^3`:

```
ell1 = 3 * omega_ob
ell2 = 3 * omega_ob^2
ell3 = omega_ob^3 * Js / npp
ell4 = 0
```

## Codebase Structure

The simulation is in `simulation/tutorials_ep6_svpwm.py`. Key components:

### 1. Controller class: `The_Motor_Controller` (jitclass)

Observer-related fields:
```python
('xS', float64[:]),        # observer states, length NS=6 (pre-allocated)
('NS', int32),             # = 6 (max state count)
('ell1', float64),         # observer gain 1
('ell2', float64),         # observer gain 2
('ell3', float64),         # observer gain 3
('ell4', float64),         # observer gain 4
('index_separate_speed_estimation', int32),      # 0=direct, 1=observer
('use_disturbance_feedforward_rejection', int32), # 0=off, 1=xS[2], 2=xS[2]+correction
('total_disrubance_feedforward', float64),        # the feedforward value to use
('speed_observer_output_error', float64),
```

Observer gains are initialized in `__init__` around line 165-185. Currently the **3rd-order** is active (`elif True`):

```python
omega_ob = 100  # [rad/s]
self.ell1 = 0.0
self.ell2 = 0.0
self.ell3 = 0.0
self.ell4 = 0.0
if False:    # 2nd-order speed observer
    self.ell2 = 2 * omega_ob
    self.ell3 = omega_ob**2 * init_Js/init_npp
elif False:  # 2nd-order position observer
    self.ell1 = 2 * omega_ob
    self.ell2 = omega_ob**2 * init_Js/init_npp
elif True:   # 3rd-order position observer  <--- CURRENTLY ACTIVE
    self.ell1 = 3 * omega_ob
    self.ell2 = 3 * omega_ob**2
    self.ell3 = omega_ob**3 * init_Js/init_npp
else:        # 4th-order position observer
    self.ell1 = 4 * omega_ob
    self.ell2 = 6 * omega_ob**2
    self.ell3 = 4 * omega_ob**3 * init_Js/init_npp
    self.ell4 = omega_ob**4
```

### 2. Observer dynamics: `DYNAMICS_SpeedObserver(x, CTRL)`

```python
@njit(nogil=True)
def DYNAMICS_SpeedObserver(x, CTRL):
    fx = np.zeros(6)
    output_error = angle_diff(CTRL.theta_d, x[0])
    CTRL.speed_observer_output_error = output_error

    fx[0] = CTRL.ell1*output_error + x[1]
    fx[1] = CTRL.ell2*output_error + (CTRL.Tem + x[2]) * CTRL.npp/CTRL.Js
    fx[2] = CTRL.ell3*output_error + x[3]
    fx[3] = CTRL.ell4*output_error + 0.0
    return fx
```

This already supports 4th-order (x[3] = p_to is driven by ell4). The dynamics are correct. The issue is only that: (a) 3rd-order gains are selected, (b) feedforward is not wired in.

### 3. Speed estimation in DSP (line 740-761)

```python
if CTRL.index_separate_speed_estimation == 0:
    CTRL.omega_r_elec = ACM.omega_r_elec          # direct measurement
elif CTRL.index_separate_speed_estimation == 1:
    RK4_ObserverSolver_CJH_Style(DYNAMICS_SpeedObserver, CTRL.xS, CTRL.CL_TS, CTRL)
    # ... angle wrapping ...
    CTRL.omega_r_elec = CTRL.xS[1]
    if CTRL.use_disturbance_feedforward_rejection == 0:
        CTRL.total_disrubance_feedforward = 0.0
    if CTRL.use_disturbance_feedforward_rejection == 1:
        CTRL.total_disrubance_feedforward = CTRL.xS[2]
    elif CTRL.use_disturbance_feedforward_rejection == 2:
        CTRL.total_disrubance_feedforward = CTRL.xS[2] + CTRL.ell2*CTRL.speed_observer_output_error
```

### 4. FOC speed controller (line 628-639)

```python
def FOC(CTRL, reg_speed, reg_id, reg_iq):
    reg_speed.setpoint = CTRL.cmd_rpm / 60 * 2*np.pi * CTRL.npp
    reg_speed.measurement = CTRL.omega_r_elec
    # ... velocity loop execution ...
    tustin_pid(reg_speed)

    if CTRL.bool_apply_speed_closed_loop_control:
        CTRL.cmd_idq[1] = reg_speed.Out    # <-- feedforward NOT added here!
```

**Problem**: `total_disrubance_feedforward` is computed but NEVER used. There's a commented-out line (685) showing the intended usage:
```python
# CTRL.cmd_idq[1] = HUMAN.KP*(reg_speed.setpoint-reg_speed.measurement) + OB.total_disrubance_feedforward
```

### 5. Load torque application

In the `user_system_input_code` string (executed via `exec()` in `start_simulation_slices`), you can set `ACM.TLoad` per time slice. The `d` dict has:
```python
'user_system_input_code': '''if ii < 1: CTRL.cmd_rpm = 50
elif ii < 5: ACM.TLoad = 0.2
elif ii < 100: CTRL.cmd_rpm = -50'''
```

### 6. Watch data for plotting

`CTRL.xS[2]` (estimated disturbance) is saved at watch index 17:
```python
watch_data[17][watch_index] = CTRL.xS[2]   # estimated TL
watch_data[18][watch_index] = CTRL.xS[3]   # estimated pT (rate of TL)
watch_data[41][watch_index] = ACM.TLoad     # actual load torque
```

And in `Watch_Mapping`:
```python
'[Nm]=CTRL.xS[2]',   # index 17: estimated -TL
'[Nm/s]=CTRL.xS[3]', # index 18: estimated DL (rate)
'[Nm]=ACM.TLoad',     # index 41: actual load
```

### 7. Simulation_Benchmark class (line 1307)

Entry point. Takes a `d` dict, auto-tunes PI if `CL_SERIES_KP` is None, runs simulation in slices.

```python
sim = Simulation_Benchmark(d)
# Access results:
sim.gdd['CTRL.xS[2]']      # estimated disturbance
sim.gdd['ACM.TLoad']        # actual load torque
sim.gdd['CTRL.omega_r_mech'] # speed
sim.global_machine_times     # time array
```

## What You Need To Do

Create a **new standalone Python script** `simulation/demo_4th_order_ESO.py` that:

### Step 1: Copy the base config

Use the `d_user_input` dict from `tutorials_ep6_svpwm.py` (line 1475) as a starting point. Use these motor parameters:
```python
'init_npp': 22,
'init_IN': 1.3*6/1.414,
'init_R': 0.035,
'init_Ld': 0.036e-3,
'init_Lq': 0.036e-3,
'init_KE': 0.0125,
'init_Rreq': 0.0,
'init_Js': 0.44e-4,
'DC_BUS_VOLTAGE': 5,
```

### Step 2: Modify the code to support observer order selection and feedforward wiring

You have two options:

**Option A (minimal, monkey-patch):** After creating `Simulation_Benchmark` with `bool_start_simulation=False`, call `get_global_objects()`, then modify `CTRL.ell1-ell4` to 4th-order gains, set `CTRL.index_separate_speed_estimation = 1`, set `CTRL.use_disturbance_feedforward_rejection = 1`, and run the simulation loop manually.

**Option B (modify source):** Edit `tutorials_ep6_svpwm.py` to:
- Add a `d` dict key like `'observer_order': 4` and use it in `get_global_objects()` to select gain formula
- Wire the feedforward into FOC: change line 639 from `CTRL.cmd_idq[1] = reg_speed.Out` to `CTRL.cmd_idq[1] = reg_speed.Out + CTRL.total_disrubance_feedforward / (1.5 * CTRL.npp * CTRL.KA)` (convert disturbance torque to iq current: `iq_ff = T_ff / (1.5 * npp * KA)`)

**Use Option B** — it is cleaner.

### Step 3: Apply ramp load torque

In `user_system_input_code`, create a ramp load:
```python
'user_system_input_code':
    "CTRL.cmd_rpm = 50\n"
    "if ii >= 2: ACM.TLoad = 0.05 * (ii - 2) * d['TIME_SLICE']"
```

This creates a linearly increasing load starting from time slice 2. Adjust the slope so the motor doesn't stall.

### Step 4: Run 4 experiments and plot

Run these 4 scenarios:
1. **3rd-order observer, feedforward OFF** — baseline
2. **3rd-order observer, feedforward ON** — can track step but not ramp
3. **4th-order observer, feedforward OFF** — better estimation, no compensation
4. **4th-order observer, feedforward ON** — should track ramp load

For each, plot:
- **Subplot 1**: Speed command vs actual speed (rpm)
- **Subplot 2**: Actual load torque `ACM.TLoad` vs estimated disturbance `CTRL.xS[2]`
- **Subplot 3**: `CTRL.xS[3]` (estimated rate of disturbance) — should converge to the ramp slope for 4th-order

Arrange as a 2x2 grid of multi-subplot figures, or a single figure with 4 columns.

### Step 5: Also test step load for comparison

Add a second experiment set with step load (`ACM.TLoad = 0.2` at some time) to show both observers can handle step, but only 4th-order handles ramp.

## Important Implementation Details

1. **Feedforward unit conversion**: The observer estimates disturbance in mechanical acceleration units (from the motion equation `d(omega)/dt = (Tem + d_to) * npp/Js`). Actually, looking at the dynamics more carefully:
   - `x[2]` has units such that `(CTRL.Tem + x[2]) * CTRL.npp / CTRL.Js` gives `d(omega_elec)/dt`
   - So `x[2]` is in the same units as `Tem` (Nm-like, but scaled)
   - The feedforward should be: `iq_ff = CTRL.total_disrubance_feedforward / (1.5 * npp * KA)`
   - But be careful: check the sign. The disturbance is `d_to = -(TLoad + B*omega)/npp_scale`, so the feedforward should SUBTRACT it. Look at the motion equation in the observer: `Tem + x[2]` means x[2] already has the right sign convention for adding to torque.

2. **Observer bandwidth `omega_ob`**: The default is 100 rad/s. For 4th-order, you may need to reduce it (e.g., 80 rad/s) because higher-order observers are more sensitive to noise. Or keep 100 and see what happens.

3. **The `numba` `@jitclass` constraint**: You CANNOT add new fields to `The_Motor_Controller` without also adding them to the `spec` list. The existing `ell1-ell4` and `xS[0:5]` already support 4th-order — no spec changes needed.

4. **The key modification to FOC()**: Around line 638-639, change:
```python
if CTRL.bool_apply_speed_closed_loop_control == True:
    CTRL.cmd_idq[1] = reg_speed.Out
```
to:
```python
if CTRL.bool_apply_speed_closed_loop_control == True:
    CTRL.cmd_idq[1] = reg_speed.Out
    if CTRL.use_disturbance_feedforward_rejection > 0:
        # Convert estimated disturbance torque to iq current feedforward
        if CTRL.KA > 0:
            CTRL.cmd_idq[1] += CTRL.total_disrubance_feedforward / (1.5 * CTRL.npp * CTRL.KA)
```

5. **Sign convention check**: In `DYNAMICS_MACHINE`, the motion equation is:
```python
fx[1] = (ACM.Tem - ACM.TLoad) / ACM.Js  # mech. angular rotor speed
```
   But in the observer:
```python
fx[1] = CTRL.ell2*output_error + (CTRL.Tem + x[2]) * CTRL.npp/CTRL.Js
```
   So `x[2]` estimates `-TLoad` (negative of load torque) scaled by npp. The feedforward `total_disrubance_feedforward = xS[2]` is negative when there's positive load. So adding it to iq command (which produces positive torque) should cancel the load. Verify this by checking the sign in your plots.

## File Dependencies

```
simulation/
  tutorials_ep6_svpwm.py   # main simulation (modify this)
  tuner.py                  # PI auto-tuning (no changes needed)
  demo_4th_order_ESO.py     # NEW FILE: your comparison script
```

## Expected Output

A matplotlib figure saved as `simulation/fig_4th_order_ESO_comparison.png` showing that:
- Under ramp load, 3rd-order observer's disturbance estimate lags behind (steady-state error in estimation)
- Under ramp load, 4th-order observer's disturbance estimate tracks the ramp with zero steady-state error
- With feedforward ON + 4th-order observer, speed tracking under ramp load is significantly better
- Under step load, both observers perform comparably (4th-order may be slightly better in transient)
