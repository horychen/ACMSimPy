#!/usr/bin/env python
"""
Staged tuning evaluation for sensorless PMSM control.

Runs Stage 0 → 3 sequentially on a given motor, collecting progressively
harder metrics.  Results are printed as a consolidated scorecard.

Usage:
    python eval_staged_tuning.py --motor servo
    python eval_staged_tuning.py --motor small_L
"""
import sys, os, copy, time, argparse
import numpy as np

# ── import the demo module ──────────────────────────────────────────
from demo_sensorless_active_flux import (
    run_sensorless_demo, plot_sensorless_results, angle_error
)

# ====================================================================
#  Metric helpers
# ====================================================================
def compute_metrics(results, cmd_rpm_ref):
    """Compute a comprehensive metric dict from simulation results."""
    t   = results['t']
    n   = len(t)
    npp = results['n_pp']

    # Define time windows -----------------------------------------------
    # Ramp-up:    0 ~ 0.3 s
    # Steady 1:   0.4 ~ 0.55 s  (before load)
    # Load step:  0.6 ~ 0.75 s  (during load)
    # Load off:   0.85 ~ 0.95 s (after load removed)
    # Reversal:   1.05 ~ 1.25 s (speed reversal transient)
    # Steady 2:   1.5 ~ 1.9 s   (reverse steady state)
    windows = {
        'ss1':      (0.40, 0.55),
        'load_on':  (0.60, 0.75),
        'load_off': (0.85, 0.95),
        'reversal': (1.05, 1.25),
        'ss2':      (1.50, 1.90),
    }

    def idx_range(t0, t1):
        return (t >= t0) & (t <= t1)

    # ── Angle metrics ─────────────────────────────────────────────────
    theta_true = results['theta_true']
    theta_af   = results['theta_af']
    err_af_deg = np.abs(angle_error(theta_af, theta_true)) * 180 / np.pi

    # ── Speed metrics ─────────────────────────────────────────────────
    omega_true = results['omega_true']  # rpm
    omega_obs  = results['omega_af_eso']  # observer speed (ESO or NSO)
    omega_lpf  = results['omega_af_lpf']
    cmd_rpm    = results['cmd_rpm']
    spd_err_obs = omega_obs - omega_true
    spd_err_lpf = omega_lpf - omega_true
    trk_err     = omega_true - cmd_rpm  # tracking error (true vs cmd)

    # ── Load torque metrics ───────────────────────────────────────────
    load_true = results['load_true']
    load_est  = results['load_eso']

    m = {}

    # ----------------------------------------------------------------
    # Level 1: Basic stability (can it run without blowing up?)
    # ----------------------------------------------------------------
    m['stable'] = not (np.any(np.isnan(omega_true)) or np.any(np.abs(omega_true) > 1e5))
    m['max_speed_abs'] = np.max(np.abs(omega_true))

    # ----------------------------------------------------------------
    # Level 2: Angle estimation accuracy (AF observer quality)
    # ----------------------------------------------------------------
    mask_ss = idx_range(*windows['ss1']) | idx_range(*windows['ss2'])
    mask_all = idx_range(0.1, 1.9)  # skip first 100ms init

    m['angle_rms_ss_deg']  = np.sqrt(np.mean(err_af_deg[mask_ss]**2)) if mask_ss.any() else 999
    m['angle_max_ss_deg']  = np.max(err_af_deg[mask_ss]) if mask_ss.any() else 999
    m['angle_rms_all_deg'] = np.sqrt(np.mean(err_af_deg[mask_all]**2)) if mask_all.any() else 999

    # ----------------------------------------------------------------
    # Level 3: Speed estimation accuracy (observer quality)
    # ----------------------------------------------------------------
    m['spd_est_rms_ss']  = np.sqrt(np.mean(spd_err_obs[mask_ss]**2)) if mask_ss.any() else 999
    m['spd_est_rms_all'] = np.sqrt(np.mean(spd_err_obs[mask_all]**2)) if mask_all.any() else 999
    m['spd_est_max_ss']  = np.max(np.abs(spd_err_obs[mask_ss])) if mask_ss.any() else 999

    # ----------------------------------------------------------------
    # Level 4: Speed tracking accuracy (controller + observer closed loop)
    # ----------------------------------------------------------------
    m['trk_rms_ss']  = np.sqrt(np.mean(trk_err[mask_ss]**2)) if mask_ss.any() else 999
    m['trk_max_ss']  = np.max(np.abs(trk_err[mask_ss])) if mask_ss.any() else 999

    # ----------------------------------------------------------------
    # Level 5: Load torque estimation
    # ----------------------------------------------------------------
    mask_load = idx_range(*windows['load_on'])
    if mask_load.any():
        m['tl_mean_true'] = np.mean(load_true[mask_load])
        m['tl_mean_est']  = np.mean(load_est[mask_load])
        m['tl_ratio']     = abs(m['tl_mean_est'] / m['tl_mean_true']) if abs(m['tl_mean_true']) > 1e-6 else 0
    else:
        m['tl_mean_true'] = m['tl_mean_est'] = m['tl_ratio'] = 0

    # ----------------------------------------------------------------
    # Level 6: Load disturbance rejection (speed dip under load step)
    # ----------------------------------------------------------------
    mask_pre  = idx_range(0.55, 0.60)
    mask_post = idx_range(0.60, 0.80)
    if mask_pre.any() and mask_post.any():
        spd_before = np.mean(omega_true[mask_pre])
        spd_worst  = omega_true[mask_post]
        m['load_dip_rpm'] = abs(spd_before - np.min(spd_worst)) if spd_before > 0 else abs(spd_before - np.max(spd_worst))
    else:
        m['load_dip_rpm'] = 999

    # ----------------------------------------------------------------
    # Level 7: Reversal transient (settling time after speed reversal)
    # ----------------------------------------------------------------
    mask_rev = idx_range(1.0, 1.9)
    if mask_rev.any():
        rev_idx = np.where(mask_rev)[0]
        target_rpm = -cmd_rpm_ref
        settled = np.abs(omega_true[rev_idx] - target_rpm) < 0.05 * cmd_rpm_ref
        if settled.any():
            first_settled = np.argmax(settled)
            m['reversal_settling_s'] = t[rev_idx[first_settled]] - 1.0
        else:
            m['reversal_settling_s'] = 999
    else:
        m['reversal_settling_s'] = 999

    return m


# ====================================================================
#  Print scorecard
# ====================================================================
def print_scorecard(stage_results, cmd_rpm_ref):
    """Print a consolidated comparison table."""
    print('\n' + '=' * 100)
    print('  STAGED TUNING SCORECARD')
    print('=' * 100)

    levels = [
        ('L1: Stable?',                       'stable',              '{:}',    None),
        ('L2: Angle RMS (SS) [deg]',          'angle_rms_ss_deg',   '{:>8.2f}', 2.0),
        ('L2: Angle RMS (all) [deg]',         'angle_rms_all_deg',  '{:>8.2f}', 5.0),
        ('L2: Angle Max (SS) [deg]',          'angle_max_ss_deg',   '{:>8.2f}', 5.0),
        ('L3: Speed est. RMS (SS) [rpm]',     'spd_est_rms_ss',    '{:>8.2f}', 0.05 * cmd_rpm_ref),
        ('L3: Speed est. RMS (all) [rpm]',    'spd_est_rms_all',   '{:>8.2f}', 0.10 * cmd_rpm_ref),
        ('L3: Speed est. Max (SS) [rpm]',     'spd_est_max_ss',    '{:>8.2f}', 0.10 * cmd_rpm_ref),
        ('L4: Tracking RMS (SS) [rpm]',       'trk_rms_ss',        '{:>8.2f}', 0.02 * cmd_rpm_ref),
        ('L4: Tracking Max (SS) [rpm]',       'trk_max_ss',        '{:>8.2f}', 0.05 * cmd_rpm_ref),
        ('L5: TL ratio (est/true)',            'tl_ratio',          '{:>8.3f}', None),
        ('L6: Load speed dip [rpm]',           'load_dip_rpm',      '{:>8.2f}', 0.10 * cmd_rpm_ref),
        ('L7: Reversal settling [s]',          'reversal_settling_s','{:>8.3f}', 0.5),
    ]

    stage_names = [s[0] for s in stage_results]
    header = f'  {"Metric":<38s}' + ''.join(f'{name:>16s}' for name in stage_names)
    print(header)
    print('  ' + '-' * (38 + 16 * len(stage_names)))

    for label, key, fmt, threshold in levels:
        row = f'  {label:<38s}'
        for _, m in stage_results:
            val = m.get(key, '—')
            if isinstance(val, bool):
                cell = 'YES' if val else 'NO '
                cell = f'{cell:>8s}'
            elif val == 999:
                cell = f'{"FAIL":>8s}'
            else:
                cell = fmt.format(val)
            # Color coding (pass/fail indicator)
            if threshold is not None and not isinstance(val, bool) and val != 999:
                marker = ' <=' if val <= threshold else ' !!'
            else:
                marker = '   '
            row += f'{cell}{marker}   '
        print(row)

    print('=' * 100)
    print(f'  Thresholds: L2 angle <2/5 deg, L3 speed <{0.05*cmd_rpm_ref:.0f}/{0.10*cmd_rpm_ref:.0f} rpm, '
          f'L4 tracking <{0.02*cmd_rpm_ref:.0f}/{0.05*cmd_rpm_ref:.0f} rpm, L7 settling <0.5 s')


# ====================================================================
#  Motor presets (mirror of main script)
# ====================================================================
def get_motor_preset(motor_name):
    common_ctrl = {
        'CTRL.bool_apply_speed_closed_loop_control': True,
        'CTRL.bool_apply_decoupling_voltages_to_current_regulation': False,
        'CTRL.bool_apply_sweeping_frequency_excitation': False,
        'CTRL.bool_overwrite_speed_commands': True,
        'CTRL.bool_zero_id_control': True,
        'FOC_CL_KI_factor_when__bool_apply_decoupling_voltages_to_current_regulation__is_False': 10,
        'CL_SERIES_KP': None, 'CL_SERIES_KI': None,
        'VL_SERIES_KP': None, 'VL_SERIES_KI': None,
        'VL_LIMIT_OVERLOAD_FACTOR': 10.0,  # generous: avoid current saturation during algorithm validation
        'disp.Kp': 0.0, 'disp.Ki': 0.0, 'disp.Kd': 0.0,
        'disp.tau': 0.0, 'disp.OutLimit': 0.0, 'disp.IntLimit': 0.0,
    }
    presets = {
        'small_L': {
            'd': {'CL_TS': 1e-4, 'VL_EXE_PER_CL_EXE': 5,
                  'MACHINE_SIMULATIONs_PER_SAMPLING_PERIOD': 1,
                  'TIME_SLICE': 0.1, 'NUMBER_OF_SLICES': 20,
                  'init_npp': 22, 'init_IN': 1.3*6/1.414,
                  'init_R': 0.035, 'init_Ld': 0.036e-3, 'init_Lq': 0.036e-3,
                  'init_KE': 0.0125, 'init_Rreq': 0.0, 'init_Js': 0.44e-4,
                  'DC_BUS_VOLTAGE': 20,  # no voltage saturation at 200 rpm
                  'FOC_delta': 15, 'FOC_desired_VLBW_HZ': 120,
                  **common_ctrl},
            'clbw': 320, 'af_kp': 500, 'af_ki': 5000, 'zeta': 15,
            'r_mis': 1.5, 'v_off': (0.02, -0.015),
            'omega_ob': 200.0,
            'cmd_rpm': 200, 'load_step': 0.15,
        },
        'servo': {
            'd': {'CL_TS': 1e-4, 'VL_EXE_PER_CL_EXE': 5,
                  'MACHINE_SIMULATIONs_PER_SAMPLING_PERIOD': 1,
                  'TIME_SLICE': 0.1, 'NUMBER_OF_SLICES': 20,
                  'init_npp': 4, 'init_IN': 4.0,
                  'init_R': 1.1, 'init_Ld': 5e-3, 'init_Lq': 6e-3,
                  'init_KE': 0.1, 'init_Rreq': 0.0, 'init_Js': 0.008,
                  'DC_BUS_VOLTAGE': 200,
                  'FOC_delta': 10, 'FOC_desired_VLBW_HZ': 60,
                  **common_ctrl},
            'clbw': 500, 'af_kp': 500, 'af_ki': 1000, 'zeta': 10,
            'r_mis': 1.5, 'v_off': (0.5, -0.3),
            'omega_ob': 80.0,
            'cmd_rpm': 500, 'load_step': 1.27,
        },
        'big_L': {
            'd': {'CL_TS': 1e-4, 'VL_EXE_PER_CL_EXE': 5,
                  'MACHINE_SIMULATIONs_PER_SAMPLING_PERIOD': 1,
                  'TIME_SLICE': 0.1, 'NUMBER_OF_SLICES': 20,
                  'init_npp': 24, 'init_IN': 4.93,
                  'init_R': 1.97, 'init_Ld': 0.1035, 'init_Lq': 0.1063,
                  'init_KE': 0.0745, 'init_Rreq': 0.0, 'init_Js': 1.5*0.051,
                  'DC_BUS_VOLTAGE': 800,
                  'FOC_delta': 6.5, 'FOC_desired_VLBW_HZ': 40,
                  **common_ctrl},
            'clbw': 200, 'af_kp': 100, 'af_ki': 1000, 'zeta': 6.5,
            'r_mis': 1.5, 'v_off': (1.0, -0.8),
            'omega_ob': 200.0,
            'cmd_rpm': 150, 'load_step': 20.0,
        },
    }
    return presets[motor_name]


# ====================================================================
#  Run one stage
# ====================================================================
def run_stage(motor_name, preset, stage, observer='nso', pure_p=False, verbose=True):
    """
    stage: 0 = encoder baseline
           1 = AF tuning (encoder control, observer passively logging)
           2 = torque feedforward only (encoder speed + observer TL)
           2b= observer speed + encoder angle (no TL FF)
           3 = full sensorless (observer speed + observer angle + TL FF)
    """
    stage_str = str(stage)
    use_speed = stage_str in ('2b', '3')
    use_angle = (stage_str == '3')  # only Stage 3 uses AF angle for Park
    use_tqff  = stage_str in ('2', '3')
    
    tag = f'Stage{stage_str}'
    if verbose:
        print(f'\n{"#"*70}')
        print(f'  {tag}: motor={motor_name}, observer={observer}, '
              f'sensorless_speed={use_speed}, sensorless_angle={use_angle}, '
              f'tqff={use_tqff}, pure_p={pure_p}')
        print(f'{"#"*70}')

    p = preset
    results = run_sensorless_demo(
        p['d'],
        zeta=p['zeta'],
        CLBW_Hz=p['clbw'],
        af_Kp=p['af_kp'],
        af_Ki=p['af_ki'],
        R_mismatch_factor=p['r_mis'],
        voltage_offset_alpha=p['v_off'][0],
        voltage_offset_beta=p['v_off'][1],
        eso_omega_ob=p.get('omega_ob', 200.0),
        speed_observer=observer,
        use_sensorless_speed=use_speed,
        use_sensorless_angle=use_angle,
        use_sensorless_torque_ff=use_tqff,
        pure_p_current=pure_p,
        cmd_rpm_ref=p['cmd_rpm'],
        load_step=p['load_step'],
        verbose=verbose,
    )

    # Save plots
    suffix = f'_s{stage}_{observer}'
    base = f'fig_eval_{motor_name}{suffix}'
    plot_sensorless_results(results, save_path=base)
    import matplotlib.pyplot as plt
    plt.close('all')

    metrics = compute_metrics(results, p['cmd_rpm'])
    return tag, metrics


# ====================================================================
#  Main
# ====================================================================
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Staged tuning evaluation')
    parser.add_argument('--motor', type=str, choices=['small_L', 'servo', 'big_L'], default='servo')
    parser.add_argument('--observer', type=str, choices=['eso', 'nso'], default='nso')
    parser.add_argument('--pure-p', action='store_true', default=False,
                        help='Use pure-P current control (stage 2+3 only)')
    parser.add_argument('--stages', type=str, default='0,1,2,2b,3',
                        help='Comma-separated stage numbers to run (default: 0,1,2,2b,3)')
    args = parser.parse_args()

    preset = get_motor_preset(args.motor)
    stages = [s.strip() for s in args.stages.split(',')]  # keep as strings for '2b'

    all_results = []
    t_start = time.time()

    for stage in stages:
        use_pure_p = args.pure_p and stage in ('2', '2b', '3')
        tag, metrics = run_stage(args.motor, preset, stage, observer=args.observer, pure_p=use_pure_p)
        all_results.append((tag, metrics))

    elapsed = time.time() - t_start
    print(f'\n  Total evaluation time: {elapsed:.1f} s')

    print_scorecard(all_results, preset['cmd_rpm'])
