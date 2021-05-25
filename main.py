# %%
from acmsimpy import *
from datetime import time
ACMSimPy(1, 0.1); # Compile

# %%
%time watch = ACMSimPy(omega_ob=1, TIME=5.5, MACHINE_TS=1e-4, CL_TS=1e-4)
times, id, iq, ia, ib, speed = watch
plt.figure(figsize=(12, 6))
plt.plot(times, id)
plt.plot(times, iq)
plt.figure(figsize=(12, 6))
plt.plot(times, ia)
plt.plot(times, ib)
plt.figure(figsize=(12, 6))
plt.plot(times, speed)
# %time (*ACMSimPy(omega_ob=1, TIME=100, MACHINE_TS=0.0001))

# %%


