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

    # Filtered derivative on measurement: negate the measurement increment only.
    reg.differentiator = ((2.0 * reg.tau - reg.T) * reg.differentiator \
                        - 2.0 * reg.Kd * (reg.measurement - reg.prevMeasurement)) \
                        / (2.0 * reg.tau + reg.T)

    # Clamp the integral to the remaining output range in this sample.
    # Separate bounds also work when P + D exceeds either output limit.
    integral_min = -reg.OutLimit - proportional - reg.differentiator
    integral_max =  reg.OutLimit - proportional - reg.differentiator
    if reg.integrator > integral_max:
        reg.integrator = integral_max
    elif reg.integrator < integral_min:
        reg.integrator = integral_min

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
