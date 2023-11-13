import numpy as np
def humans_give_commands(CTRL, ACM, t):
    """ Console @ CL_TS """

    if t < 1:
        CTRL.cmd_rpm = 100
    elif t < 2:
        CTRL.cmd_rpm = 200
    # elif t < 0.4:
    #     CTRL.cmd_rpm = 50
    # elif t < 0.6:
    #     CTRL.cmd_rpm = 100
    # elif t < 0.8:
    #     CTRL.cmd_rpm = 150
    # elif t < 1.0:
    #     CTRL.cmd_rpm = 200
    # elif t < 1.2:
    #     CTRL.cmd_rpm = 250
    # elif t < 1.4:
    #     CTRL.cmd_rpm = 700
    # elif t < 1.6:
    #     CTRL.cmd_rpm = 1200
    # elif t < 1.8:
    #     CTRL.cmd_rpm = 1500
    elif t < 3:
         ACM.TLoad = 2
    # elif t < 2.0:
    #      ACM.TLoad = 4
    # elif t < 6: 
    #      CTRL.CMD_SPEED_SINE_RPM = 200s
    # elif t < 2.1:
    #      ACM.TLoad = 4
    # elif t < 2.2:
    #      ACM.TLoad = 4.77
    # elif t < 2.3:
    #     CTRL.cmd_rpm = 500
    # elif t < 4:
    #     CTRL.cmd_rpm = 600
    elif t < 4:
        CTRL.cmd_rpm = 300
    elif t < 2:
        
        CTRL.apply_pulse_4_evaluating_position_estimator_accuracy = False
    # elif t < 2.40:
    #     CTRL.apply_pulse_4_evaluating_position_estimator_accuracy = False
    # elif t < 4:
    #     # CTRL.cmd_rpm = -200
    #     CTRL.apply_pulse_4_evaluating_position_estimator_accuracy = True
    # elif t < 4.1:
    #     # CTRL.cmd_rpm = -200
    #     CTRL.apply_pulse_4_evaluating_position_estimator_accuracy = True
    # elif t < 4.0:
    # print(CTRL.cmd_rpm)

    if CTRL.bool_overwrite_speed_commands == False:
        if t < 1.0:
            CTRL.cmd_rpm = 50
        elif t < 1.5:
            ACM.TLoad = 2
        elif t < 2.0:
            CTRL.cmd_rpm = 200
        elif t < 3.0:
            CTRL.cmd_rpm = -200
        elif t < 4.0:
            CTRL.cmd_rpm = 0
        elif t < 4.5:
            CTRL.cmd_rpm = 2000
        elif t < 5:
            CTRL.cmd_idq[0] = 2
        elif t < 5.5:
            ACM.TLoad = 0.0
        elif t < 6: 
            CTRL.CMD_SPEED_SINE_RPM = 500
        # else: # don't implement else to receive commands from IPython console

        # if CTRL.CMD_SPEED_SINE_RPM!=0:
        #     CTRL.cmd_rpm = CTRL.CMD_SPEED_SINE_RPM * np.sin(2*np.pi*CTRL.CMD_SPEED_SINE_HZ*t)
        pass

    if CTRL.bool_apply_sweeping_frequency_excitation == True:

        if CTRL.timebase > CTRL.CMD_SPEED_SINE_END_TIME:
            # next frequency
            CTRL.CMD_SPEED_SINE_HZ += CTRL.CMD_SPEED_SINE_STEP_SIZE
            # next end time
            CTRL.CMD_SPEED_SINE_LAST_END_TIME = CTRL.CMD_SPEED_SINE_END_TIME
            CTRL.CMD_SPEED_SINE_END_TIME += 1.0/CTRL.CMD_SPEED_SINE_HZ # 1.0 Duration for each frequency

        if CTRL.CMD_SPEED_SINE_HZ > CTRL.CMD_SPEED_SINE_HZ_CEILING:
            # stop
            CTRL.cmd_rpm = 0.0
            CTRL.cmd_idq[1] = 0.0
        else:
            # speed control - closed-loop sweep
            CTRL.cmd_rpm    = CTRL.CMD_SPEED_SINE_RPM      * np.sin(2*np.pi*CTRL.CMD_SPEED_SINE_HZ*(CTRL.timebase - CTRL.CMD_SPEED_SINE_LAST_END_TIME))

            # speed control - open-loop sweep
            CTRL.cmd_idq[1] = CTRL.CMD_CURRENT_SINE_AMPERE * np.sin(2*np.pi*CTRL.CMD_SPEED_SINE_HZ*(CTRL.timebase - CTRL.CMD_SPEED_SINE_LAST_END_TIME))

    if CTRL.bool_yanzhengzhang == True:
        # dasdasdsa
        pass
