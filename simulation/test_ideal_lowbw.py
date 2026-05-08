#!/usr/bin/env python
"""Test with IDEAL parameters: no R mismatch, no voltage offset."""
import sys; sys.path.insert(0, '.')
from test_alt_observers import run_alt_observer_sim, eval_results
from demo_sensorless_active_flux import run_sensorless_demo
from eval_staged_tuning import get_motor_preset, compute_metrics
import copy

for motor_name in ['servo', 'big_L']:
    p = get_motor_preset(motor_name)
    cmd = p['cmd_rpm']
    R = p['d']['init_R']; L = p['d']['init_Lq']
    print(f'\n{"="*70}')
    print(f'  {motor_name} IDEAL params (R_mis=1.0, v_off=0)')
    print(f'{"="*70}')
    
    # AF with ideal params
    print('\n  --- AF (ideal) ---')
    clbws = [10, 20, 30, 50, 100] if motor_name == 'servo' else [5, 10, 20, 50]
    for clbw in clbws:
        try:
            r = run_sensorless_demo(
                p['d'], zeta=p['zeta'], CLBW_Hz=clbw,
                af_Kp=p['af_kp'], af_Ki=p['af_ki'],
                R_mismatch_factor=1.0, voltage_offset_alpha=0, voltage_offset_beta=0,
                eso_omega_ob=50, speed_observer='nso',
                use_sensorless_speed=True, use_sensorless_angle=True,
                use_sensorless_torque_ff=True,
                cmd_rpm_ref=cmd, load_step=p['load_step'], verbose=False)
            m = compute_metrics(r, cmd)
            rev = f'{m["reversal_settling_s"]:.3f}' if m['reversal_settling_s']<999 else 'FAIL'
            print(f'    CLBW={clbw:>4d} | Ang={m["angle_rms_ss_deg"]:>6.1f} Trk={m["trk_rms_ss"]:>6.1f} {rev}')
        except Exception as e:
            print(f'    CLBW={clbw:>4d} | ERROR')
    
    # PLL with ideal params
    print('\n  --- PLL (ideal) ---')
    configs = [(10,30), (10,50), (20,30), (20,50), (20,100),
               (30,50), (30,100), (50,100), (50,200)] if motor_name == 'servo' else \
              [(5,10), (5,20), (8,20), (10,20), (10,50), (20,50)]
    for clbw, pll_bw in configs:
        try:
            r = run_alt_observer_sim(
                p['d'], observer_type='pll', CLBW_Hz=clbw, zeta=p['zeta'],
                R_mismatch_factor=1.0, v_off=(0,0),
                omega_ob=50, cmd_rpm=cmd, load_step=p['load_step'],
                pll_bw=pll_bw, af_Kp=p['af_kp'], af_Ki=p['af_ki'], verbose=False)
            ang, spd, trk, dip, rev = eval_results(r, cmd)
            rev_s = f'{rev:.3f}' if rev<999 else 'FAIL'
            print(f'    CL={clbw:>3d} pll={pll_bw:>3d} | Ang={ang:>6.1f} Trk={trk:>6.1f} {rev_s}')
        except Exception as e:
            print(f'    CL={clbw:>3d} pll={pll_bw:>3d} | ERROR')
    
    # SMO with ideal params
    print('\n  --- SMO (ideal) ---')
    configs_smo = [(10,30), (20,50), (30,50), (30,100), (50,100)] if motor_name == 'servo' else \
                  [(5,10), (8,10), (8,20), (10,20), (10,50)]
    for clbw, pll_bw in configs_smo:
        try:
            r = run_alt_observer_sim(
                p['d'], observer_type='smo', CLBW_Hz=clbw, zeta=p['zeta'],
                R_mismatch_factor=1.0, v_off=(0,0),
                omega_ob=50, cmd_rpm=cmd, load_step=p['load_step'],
                pll_bw=pll_bw, smo_lpf=max(100, clbw*3),
                af_Kp=p['af_kp'], af_Ki=p['af_ki'], verbose=False)
            ang, spd, trk, dip, rev = eval_results(r, cmd)
            rev_s = f'{rev:.3f}' if rev<999 else 'FAIL'
            print(f'    CL={clbw:>3d} pll={pll_bw:>3d} | Ang={ang:>6.1f} Trk={trk:>6.1f} {rev_s}')
        except Exception as e:
            print(f'    CL={clbw:>3d} pll={pll_bw:>3d} | ERROR')
