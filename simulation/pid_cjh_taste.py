def incremental_pi(reg):
    reg.Err = reg.setpoint - reg.measurement
    reg.Out = reg.OutPrev + \
        reg.Kp * (reg.Err - reg.ErrPrev) + \
        reg.Ki * reg.Err
    if reg.Out >    reg.OutLimit:
        reg.Out =   reg.OutLimit
    elif reg.Out < -reg.OutLimit:
        reg.Out =  -reg.OutLimit
    reg.ErrPrev = reg.Err
    reg.OutPrev = reg.Out

def tustin_pid(reg):

    # Error signal
    error = reg.setpoint - reg.measurement

    # Proportional
    proportional = reg.Kp * error

    # Filtered derivative on measurement: negate the measurement increment only.
    reg.differentiator = ((2.0 * reg.tau - reg.T) * reg.differentiator \
                        - 2.0 * reg.Kd * (reg.measurement - reg.prevMeasurement)) \
                        / (2.0 * reg.tau + reg.T)

    # Ki = 0 disables integral action, including previously stored state.
    if reg.Ki == 0.0:
        reg.integrator = 0.0
    else:
        integral_step = 0.5 * reg.Ki * reg.T * (error + reg.prevError) # Tustin
        integral_candidate = reg.integrator + integral_step
        output_candidate = proportional + integral_candidate + reg.differentiator
        # Freeze only updates that deepen saturation. Use the actual Tustin
        # increment, whose sign can differ from the current error after reversal.
        if not ((output_candidate > reg.OutLimit and integral_step > 0.0) or
                (output_candidate < -reg.OutLimit and integral_step < 0.0)):
            reg.integrator = integral_candidate

    # Compute output and apply limits
    reg.Out = proportional + reg.integrator + reg.differentiator

    if reg.Out  >  reg.OutLimit:
        reg.Out =  reg.OutLimit
    elif reg.Out< -reg.OutLimit:
        reg.Out = -reg.OutLimit

    # Store error and measurement for later use */
    reg.prevError       = error
    reg.prevMeasurement = reg.measurement

    # Return controller output */
    return reg.Out
