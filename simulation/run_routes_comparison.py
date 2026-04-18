import time as _time
from demo_routes_comparison import get_base_d, run_scenario, plot_bode

if __name__ == '__main__':
    print('=' * 70)
    print(' 探索两条减小抗扰峰值的路线')
    print('=' * 70)

    t_total_start = _time.time()
    
    omega_ob = 500.0
    cmd_rpm = 0.0
    
    sweep_results = []
    
    # Baseline
    print('--- Running Baseline ---')
    d_base = get_base_d()
    d_base['NUMBER_OF_SLICES'] = 1300
    d_base['override_speed_ki'] = None # Auto KI, not 0
    res_base = run_scenario(d_base, observer_order=4, feedforward_on=True, load_type='sweep', omega_ob=omega_ob, cmd_rpm=cmd_rpm, verbose=False)
    res_base['label'] = 'Baseline (VLBW=120Hz, Kd=0)'
    sweep_results.append(res_base)
    
    # Route 1: High VLBW
    print('--- Running Route 1: High PI Bandwidth ---')
    d_r1 = get_base_d()
    d_r1['FOC_desired_VLBW_HZ'] = 450
    d_r1['NUMBER_OF_SLICES'] = 1300
    d_r1['override_speed_ki'] = None
    res_r1 = run_scenario(d_r1, observer_order=4, feedforward_on=True, load_type='sweep', omega_ob=omega_ob, cmd_rpm=cmd_rpm, verbose=False)
    res_r1['label'] = 'Route 1: High PI Bandwidth (VLBW=450Hz)'
    sweep_results.append(res_r1)
    
    # Route 2: Active Virtual Inertia
    print('--- Running Route 2: Active Virtual Inertia ---')
    d_r2 = get_base_d()
    d_r2['NUMBER_OF_SLICES'] = 1300
    d_r2['override_speed_ki'] = None
    
    # Init_js = 0.44e-4. Virtual Inertia Kd acts as J_virtual. 
    # setting it to 0.01 significantly increases effective inertia.
    res_r2 = run_scenario(d_r2, observer_order=4, feedforward_on=True, load_type='sweep', omega_ob=omega_ob, cmd_rpm=cmd_rpm, verbose=False, virtual_inertia_kd=0.01)
    res_r2['label'] = 'Route 2: Active Inertia (VLBW=120Hz, Kd=0.01)'
    sweep_results.append(res_r2)
    
    print('Generating Bode plots...')
    plot_bode(sweep_results, omega_ob=omega_ob, save_path='fig_routes_comparison_bode.png')
    print('Done! Saved figure to fig_routes_comparison_bode.png')
