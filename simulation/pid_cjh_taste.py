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

    # Integral
    reg.integrator = reg.integrator + 0.5 * reg.Ki * reg.T * (error + reg.prevError) # Tustin
    # reg.integrator = reg.integrator + reg.Ki * reg.T * (error) # Euler

    # Anti-wind-up via integrator clamping */
    if reg.integrator  >  reg.IntLimit:
        reg.integrator =  reg.IntLimit
    elif reg.integrator< -reg.IntLimit:
        reg.integrator = -reg.IntLimit

    # Derivative (band-limited differentiator) # Note: derivative on measurement, therefore minus sign in front of equation! */
    reg.differentiator = -(2.0 * reg.Kd * (reg.measurement - reg.prevMeasurement) \
                        + (2.0 * reg.tau - reg.T) * reg.differentiator) \
                        / (2.0 * reg.tau + reg.T)

    # Compute output and apply limits
    reg.Out = proportional + reg.integrator + reg.differentiator

    if reg.Out  >  reg.OutLimit:
        reg.Out =  reg.OutLimit
    elif reg.Out< -reg.OutLimit:
        reg.Out = -reg.OutLimit

    # Store error and measurement for later use */
    reg.prevError       = error
    reg.prevMeasurement = reg.measurement

    # Implement dynamic clamping
    reg.IntLimit = reg.OutLimit - proportional 

    # Return controller output */
    return reg.Out
