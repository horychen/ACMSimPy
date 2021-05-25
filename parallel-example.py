# %%
from numba import njit, prange, jitclass, float64, int32
from pylab import np
from datetime import time

# %%
@njit(nogil=True)
def friction_fn(v, vt):
    if v > vt:
        return - v * 3
    else:
        return - vt * 3 * np.sign(v)


@jitclass(
    spec=[
        ('value', float64),
        ('list', float64[:]),
    ]
)
class Parameter_Holder:
    def __init__(self):
        self.value = 100
        self.list = np.zeros(3, dtype=np.float64)


@njit(nogil=True)
def simulate_spring_mass_funky_damper(x0, T=10, dt=0.0001, vt=1.0):
    times = np.arange(0, T, dt)
    positions = np.zeros_like(times)

    # param = Parameter_Holder()

    v = 0
    a = 0
    x = x0
    positions[0] = x0/x0

    for ii in range(len(times)):
        if ii == 0:
            continue
        t = times[ii]
        # a = friction_fn(v, vt) - param.value*x + param.list[0] + param.list[1]
        a = friction_fn(v, vt) - 100*x;
        v = v + a*dt
        x = x + v*dt
        positions[ii] = x/x0
    return times, positions


# compile
_ = simulate_spring_mass_funky_damper(0.1)

@njit(nogil=True, parallel=True)
def run_sims(end=1000):
    for x0 in prange(int(end/0.1)):
        if x0 == 0:
            continue
        simulate_spring_mass_funky_damper(x0*0.1)
run_sims()



# %%
%%time
run_sims()
run_sims.parallel_diagnostics(level=4)

# %%
