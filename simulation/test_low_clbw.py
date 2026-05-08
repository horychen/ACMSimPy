#!/usr/bin/env python
"""
Retry PLL/SMO on servo & big_L with VERY low CLBW (10-80 Hz).
Key insight from user: high control bandwidth destabilizes sensorless control.
"""
import sys; sys.path.insert(0, '.')
from test_alt_observers import run_alt_observer_sim, eval_results
from eval_staged_tuning import get_motor_preset
import copy

def run_sweep(motor_name):
    p = get_motor_preset(motor_name)
    cmd = p['cmd_rpm']
    R = p['d']['init_R']; L = p['d']['init_Lq']; KE = p['d']['init_KE']
    print(f'\n{"="*75}')
    print(f'  {motor_name}: KE/L={KE/L:.1f}, L/R={L/R*1e3:.1f}ms, R={R}, L={L*1e3:.2f}mH')
    print(f'  Plant pole f_plant = R/(2piL) = {R/(2*3.14159*L):.1f} Hz')
    print(f'{"="*75}')

    # For servo: f_plant = R/(2πL) = 1.1/(2π*5e-3) = 35 Hz
    # For big_L: f_plant = 1.97/(2π*0.1035) = 3 Hz
    # So CLBW should be near or below f_plant for stability

    configs = []
    if motor_name == 'servo':
        # CLBW range: 10-100, focus around f_plant=35 Hz
        for obs in ['pll', 'smo']:
            for clbw in [10, 15, 20, 25, 30, 40, 50, 80]:
                for pll_bw in [30, 50, 100, 200]:
                    if pll_bw > clbw * 5:  # skip unreasonable combos
                        continue
                    zeta = 10
                    w_ob = 50  # lower NSO bandwidth too
                    smo_lpf = max(100, clbw * 3) if obs == 'smo' else None
                    label = f'{obs.upper()} CL{clbw} p{pll_bw}'
                    configs.append((obs, clbw, pll_bw, None, smo_lpf, w_ob, zeta, label))
    else:  # big_L, f_plant = 3 Hz
        for obs in ['pll', 'smo']:
            for clbw in [3, 5, 8, 10, 15, 20, 30, 50]:
                for pll_bw in [10, 20, 50, 100]:
                    if pll_bw > clbw * 5:
                        continue
                    zeta = 6.5
                    w_ob = 30
                    smo_lpf = max(50, clbw * 3) if obs == 'smo' else None
                    label = f'{obs.upper()} CL{clbw} p{pll_bw}'
                    configs.append((obs, clbw, pll_bw, None, smo_lpf, w_ob, zeta, label))

    header = f'{"Label":<20s} | {"AngSS":>7s} {"SpdSS":>7s} {"TrkSS":>7s} {"Dip":>7s} {"Rev":>6s}'
    print(header); print('-' * 65)

    best_ang = 999; best_label = ''
    for obs, clbw, pll_bw, sg, sl, w_ob, z, label in configs:
        try:
            res = run_alt_observer_sim(
                p['d'], observer_type=obs, CLBW_Hz=clbw, zeta=z,
                R_mismatch_factor=p['r_mis'], v_off=p['v_off'],
                omega_ob=w_ob, cmd_rpm=cmd, load_step=p['load_step'],
                pll_bw=pll_bw, smo_gain=sg, smo_lpf=sl,
                af_Kp=p['af_kp'], af_Ki=p['af_ki'], verbose=False)
            ang, spd, trk, dip, rev = eval_results(res, cmd)
            rev_s = f'{rev:.3f}' if rev < 999 else 'FAIL'
            flag = ' ★' if ang < 10 else (' ◆' if ang < 30 else '')
            print(f'{label:<20s} | {ang:>7.2f} {spd:>7.2f} {trk:>7.2f} {dip:>7.2f} {rev_s:>6s}{flag}')
            if ang < best_ang:
                best_ang = ang; best_label = label
        except Exception as e:
            print(f'{label:<20s} | ERROR: {str(e)[:42]}')

    print(f'\n  Best: {best_label} -> Angle RMS = {best_ang:.2f} deg')
    return best_ang

# Also test AF with very low CLBW for comparison
def run_af_lowbw(motor_name):
    from demo_sensorless_active_flux import run_sensorless_demo
    from eval_staged_tuning import compute_metrics
    p = get_motor_preset(motor_name)
    cmd = p['cmd_rpm']
    
    print(f'\n  --- AF baseline with low CLBW (Stage 3) ---')
    clbws = [10, 15, 20, 30, 50] if motor_name == 'servo' else [3, 5, 8, 10, 15, 20]
    header = f'{"CLBW":>6s} | {"AngSS":>7s} {"TrkSS":>7s} {"Dip":>7s} {"Rev":>6s}'
    print(header); print('-' * 45)
    
    for clbw in clbws:
        try:
            results = run_sensorless_demo(
                p['d'], zeta=p['zeta'], CLBW_Hz=clbw,
                af_Kp=p['af_kp'], af_Ki=p['af_ki'],
                R_mismatch_factor=p['r_mis'],
                voltage_offset_alpha=p['v_off'][0], voltage_offset_beta=p['v_off'][1],
                eso_omega_ob=50, speed_observer='nso',
                use_sensorless_speed=True, use_sensorless_angle=True,
                use_sensorless_torque_ff=True,
                cmd_rpm_ref=cmd, load_step=p['load_step'], verbose=False)
            m = compute_metrics(results, cmd)
            rev_s = f'{m["reversal_settling_s"]:.3f}' if m['reversal_settling_s'] < 999 else 'FAIL'
            flag = ' ★' if m['angle_rms_ss_deg'] < 10 else ''
            print(f'{clbw:>6d} | {m["angle_rms_ss_deg"]:>7.2f} {m["trk_rms_ss"]:>7.02f} {m["load_dip_rpm"]:>7.02f} {rev_s:>6s}{flag}')
        except Exception as e:
            print(f'{clbw:>6d} | ERROR: {str(e)[:42]}')

print('='*75)
print('  LOW CLBW SWEEP: PLL, SMO, and AF on servo & big_L')
print('='*75)

for motor in ['servo', 'big_L']:
    run_af_lowbw(motor)
    run_sweep(motor)
