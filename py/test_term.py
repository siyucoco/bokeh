''' Present an interactive function explorer with slider widgets.
There are 5 input parameters that the user can play with it
Use the ``bokeh serve`` command to run the example by executing:
    bokeh serve --show dac.py
at your command prompt. If default port is taken, you can specify
port using ' --port 5010' at the end of the bokeh command.
'''
import math
from scipy.integrate import solve_ivp, odeint
from bokeh.io import save, curdoc
from bokeh.layouts import column, row
from bokeh.model import Model
from bokeh.models import CustomJS, Slider, Callback, HoverTool, Button
from bokeh.plotting import ColumnDataSource, figure, show
from bokeh.models.widgets import Panel, Tabs
import numpy as np

# --------------------- Static Parameters    --------------------- #

b0 = 93 * (10**-5)  # 93      unit : 1/Pa
deltH_0 = 95300  #               unit: j/mol 
Tw = T_in = 298.0  #  ambient temperature,  also inlet temperature, in kelvin  unit: kelvin
T0 = 353.15  # reference temeperature to be used in the Toth isotherm   unit: kelvin
t_h0 = .37  # heterogeneity constant 
apha = 0.33
chi = 0.0
q_s0 = 3.40 # qs_var = q_s0 = 3.4 due to chi = 0 mol/kg
# R = 8.314* (10**3) # Universal gas constant - LPa/molK
Rg = .0821 # Universal gas constant in l-atm/mol/K
kT0 = 3.5 * (10 ** -2)  # used to compute kT for r_CO2... in mol/Kg-pa-sec
EaCO2 = 15200 # activation for KT -- J/mol
ps = 880.0
deltH_co2 = 75000.0  # calculate temeprature change   unit: jol/mol
# T_in = 293 
R_constant = 8.314 # jol/kelvin-mol

# ------------------ For Equation : Enegergy Ballance  -------------- #
pg = 1.87  # 
h = 13.8
Cp_g = 846.0  # J/kgK
Cp_s = 1500.0  # J/kgK

def cube(x):
    if 0<=x: return x**(1./3.)
    return -(-x)**(1./3.)

def ener_balan(v0, theta, deltZ):  # replace v0  * pg* Cp_g / (theta * deltZ)
    ener_balan_part1 = v0 * pg * Cp_g
    # print(f"v0 * pg * Cp_g: ", {ener_balan_part1})
    return (ener_balan_part1 / (theta * deltZ))

def ener_balan2(episl_r):
    ener_balan2 = (1 - episl_r) * ps * deltH_co2
    # print(f"(1 - episl_r) * ps * deltH_co2: ", {ener_balan2})
    return (ener_balan2)

def ener_balan3(a_s, Tn):
    ener_balan3 = a_s * h * (Tw - Tn)
    # print(f"a_s * h * (Tw - Tn): " , {ener_balan3})
    return (ener_balan3)

# Equation 1 Mass Balance : find co2_n

def mass_balan(v0, episl_r, deltZ):
    mass_balan = v0 / (episl_r * deltZ)
    # print(f"v0 / (episl_r * deltZ): ", {mass_balan})
    return (mass_balan)

def masss_balan2(episl_r, ps):
    mass_balan2 = (1 - episl_r) * ps
    # print(f"(1 - episl_r) * ps: ", {mass_balan2})
    return (mass_balan2)

def b(T): # will be called inside deriv
    # print(f"T", {T})
    b = b0 *(math.exp(((deltH_0 / (R_constant * T0)) * (T0 / T - 1))))
    # print(f'b, {b}')
    return b

def t_h(T): # will be call inside  deriv
    # print(T)
    # print(f'T, {T}')
    t_h = t_h0 + (apha * (1 - (T0 / T)) )
    # t_h = .37 + ( .33 * (1 - (353.15 / T)) )
    # print(t_h)
    return (t_h)

# -------------- rco2 in mol/kg -sec ------------- #
def R_co2(T, c_co2, q, b_var, t_var): 
    kn = kT0 * ( math.exp(( (-EaCO2) / (R_constant*T) )))
    # print(f"T in rco2", {T})
    print(f"kn",{kn})
    rco2_1 = Rg * T * c_co2 
    print(f"Rg * T * c_co2", {rco2_1})
    # print(f"th", {t_var})
   
    # print(q)
    rco2_2 = ((1 - ((q / q_s0) ** (t_var))) ** (1 / t_var))
    # rco2_2 = ((1 - ((q / q_s0) ** (t_var))) ** (1 / t_var))
    # print(f"rco2_2", {rco2_2})abs
    rco2_3 = q / (b_var * q_s0)
    # print(f"rco2_3", {rco2_3})
    # r_co2_part1 = rco2_1    
    rco2 = kn * (rco2_1 * rco2_2 - rco2_3)
    print(f"rco2 in rco2", {rco2})
    # r_co2_part1 =(Rg * T * c_co2 * ((1 - ((q / q_s0) ** (t_var))) ** (1 / t_var)) - q / (b_var * q_s0))
    # print(f"rco2_part1",{r_co2_part1})
    # r_co2 = kn *r_co2_part1
    # print(f'r_co2, {r_co2}')
    # print(f"b_var * qs_var, {rco2_term5}.")
    return rco2
    # rco2 in 
V = 200.0  # volume
T = 293.0 # +273 ambrient temperature
c_co2_0 = .016349 # mol/m^3    
episl_r = 0.3  # void
volumetric_flow = 5 # m^3/s

r = cube(V/(20*math.pi))
v0 = volumetric_flow / (math.pi *r*r )
L = V / (math.pi * (r ** 2))
deltZ = L / 5.0  # 5 boxes in total
a_s = 2 / r
theta = (1 - episl_r) * ps * Cp_s + episl_r * pg * Cp_g
t_var = t_h(T)
    # print(t_var)
b_var = b(T)
q_init = 0
print("b_var", {b_var})
print("t_var", {t_var})
print(deltZ)
print(f"v0", {v0})
print(f"theta", {theta})


r_co2 = R_co2(Tw, c_co2_0, 0,  b_var, t_var)
temperature_constant = ((v0  * pg* Cp_g) / (theta * deltZ))
temperature_constant2 = (1 - episl_r) * ps * deltH_co2 /theta 
temperature_constant3 = a_s * h /theta

concentration_constant = v0 / (episl_r * deltZ)
concentration_constant2 = (1 - episl_r) * ps/episl_r 

print(f"concentration", concentration_constant * c_co2_0)
print(f"rco2(1-epis)", {r_co2*concentration_constant2})

# print(f"rco2",{R_co2(Tw, c_co2_0, 0,  b_var, t_var)})
# print(f"v0pgcpg", {-temperature_constant * Tw})
# print("\n")
# print(f"1-episl_R,", {temperature_constant2 * R_co2(Tw, c_co2_0, 0,  b_var, t_var)})
# print("\n")
# print(f"a_S", {temperature_constant3})
T1 = -temperature_constant* Tw + temperature_constant* Tw + temperature_constant2* (
        R_co2(293, c_co2_0, 0,  b_var, t_var))+ temperature_constant3*(Tw - Tw) # Last Tw should be T1
print(f"T1", {T1})
co2_1dot = -concentration_constant * c_co2_0 + concentration_constant * c_co2_0 - (
        R_co2(293, c_co2_0, 0,  b_var, t_var)) * concentration_constant2
print(f"co2_dot1", {co2_1dot})
q1dot = R_co2(293, c_co2_0, 0,  b_var, t_var)
print(f"q1dot", q1dot)

print(3.5*10 **3 * math.exp(-15200 / (8.314 * 298)) * (400 / (1*(10**6))))