![acmsimpy logo](https://github.com/horychen/ACMSimPy/blob/gui/gui/images/svg_images/logo_home_acmsimpy.svg?raw=true)

# Motor Control Simulation Visualized by PyOneDark Qt Modern GUI

## installation
- 1. Install _Anaconda 3_ and open command prompt cmd.exe on Windows
- 2. `conda create -n main python=3.10`
- 3. `conda activate main`
- 4. `pip install PySide6 matplotlib pandas numba qtconsole control`
- 5. `cd /d D:\acmsimpy`
- 6. `python main.py`

> To avoid PyQt plugin issues, please go with python 3.10.

## Examples

A dual-rotor-topology axial-flux in-wheel motor is driven by a three-phase voltage source inverter. Space-vector PWM considering deadtime and DC bus are as well simulated. 

![Example 1](https://github.com/horychen/ACMSimPy/blob/numba_demo_fulldynamics_svpwm/gui/images/acmsimpy-example01.png?raw=true)


## Features

- A [unified ac motor model based on active flux](https://ieeexplore.ieee.org/document/9853634/) is implemented.
- Numerical simulation is done by Runge–Kutta method (RK4, ode4).
- Space-vector PWM is implemented according to TI ControlSUITE.
- User input logic separation from code and a modern GUI.
- Batch making technical-paper-ready figures. *The code is the figure itself not the .pdf file!*
- Real time simulation like debugging in CCS with actual motor. *This is a feature that even Simulink does not provide.*
