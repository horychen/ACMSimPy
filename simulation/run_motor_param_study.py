import time as _time
import numpy as np
import copy
from demo_4th_order_ESO import get_base_d, run_scenario, plot_bode

if __name__ == '__main__':
    print('=' * 70)
    print(' 探索不同电机本体参数对四阶 ESO 频响（特别是抗扰高峰）的影响')
    print('=' * 70)

    t_total_start = _time.time()
    
    omega_ob = 500.0
    cmd_rpm = 0.0
    slices = 1300
    
    sweep_results = []
    
    # Pre-calculated auto-tuned baseline params for VLBW=120Hz
    base_vl_kp = 0.0030213409476208927
    base_vl_ki = 41.54343802978728
    base_cl_kp = 0.33650184804127703
    base_cl_ki = 972.2222222222223

    def add_scenario(label_name, param_dict, scale_dict):
        print(f'\\n--- Running: {label_name} ---')
        d_mod = get_base_d()
        d_mod['NUMBER_OF_SLICES'] = slices
        d_mod['skip_tuning'] = True
        
        # Disable auto tuner
        d_mod['override_speed_ki'] = base_vl_ki * scale_dict.get('vl_ki', 1.0)
        d_mod['VL_SERIES_KP'] = base_vl_kp * scale_dict.get('vl_kp', 1.0)
        d_mod['VL_SERIES_KI'] = d_mod['override_speed_ki']
        
        d_mod['CL_SERIES_KP'] = base_cl_kp * scale_dict.get('cl_kp', 1.0)
        d_mod['CL_SERIES_KI'] = base_cl_ki * scale_dict.get('cl_ki', 1.0)
        
        for k, v in param_dict.items():
            d_mod[k] = v
            print(f"Modifying {k} = {v}")
            
        res = run_scenario(d_mod, observer_order=4, feedforward_on=True, 
                           load_type='sweep', omega_ob=omega_ob, cmd_rpm=cmd_rpm, verbose=False)
        res['label'] = label_name
        sweep_results.append(res)
        
    # Baseline
    add_scenario('Baseline (J_s=0.44e-4)', {}, {})
    
    # Run 1: Double Inertia (J_s x 5). Huge inertia adds enormous physical rejection.
    # Kp scales with J_s to maintain the exact same speed loop physical bandwidth.
    add_scenario('Inertia x5 (J_s x5)', 
                 {'init_Js': 0.000044 * 5}, 
                 {'vl_kp': 5.0})
    
    # Run 2: Double Flux (KE x 2). Stronger magnet, higher torque per amp.
    # Kp inversely proportional to KE to maintain same speed loop bandwidth.
    add_scenario('Double Flux (K_E x2)', 
                 {'init_KE': 0.016 * 2}, 
                 {'vl_kp': 0.5})
    
    # Run 3: Half Inductance
    # CL Kp scales with L. CL Ki scales with 1/L.
    add_scenario('Half Inductance (L / 2)', 
                 {'init_Ld': 0.0025, 'init_Lq': 0.0025}, 
                 {'cl_kp': 0.5, 'cl_ki': 2.0})
    

    print('\\nGenerating Bode plots...')
    plot_bode(sweep_results, omega_ob=omega_ob, save_path='fig_motor_param_study.png')
    print('Done! Saved figure to fig_motor_param_study.png')
