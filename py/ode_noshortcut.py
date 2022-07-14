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

b0 = 93 * (10**-5)  # 93      unit : 1/bars
deltH_0 = 95300  #               unit: j/mol 
Tw = 278.0  #  room temperature  unit: kelvin
T0 = 353.15  # temeperature in   unit: kelvin
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

R_constant = 8.314 # jol/kelvin-mol

# ------------------ For Equation : Enegergy Ballance  -------------- #
pg = 1.87  # 
h = 13.8
Cp_g = 846.0  # J/kgK
Cp_s = 1500.0  # J/kgK

# ------------------ ODE Repetitive Shortcut -------------- #

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

def cube(x):
    if 0<=x: return x**(1./3.)
    return -(-x)**(1./3.)

# ------------------ Equastions to calculate Rco2 -------------- #
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

# Calculate rco2_n (not ode)
def R_co2(T, c_co2, q, b_var, t_var):

    kn = kT0 * ( math.exp(( (-EaCO2) / (R_constant*T) )))
    # print(f"T in rco2", {T})
    # print(f"kn",{kn})
    rco2_1 = Rg * T * c_co2 
    # if t_var <0:
    #     print("alert, t_var is negative")
    if q<0:
        print(f"q in rco2", {q})
        print("alert, q is nagative")
    rco2_2 = ((1 - ((q / q_s0) ** (t_var))) ** (1 / t_var))
    # rco2_2 = ((1 - ((q / q_s0) ** (t_var))) ** (1 / t_var))
    # print(f"rco2_2", {rco2_2})abs
    rco2_3 = q / (b_var * q_s0)
    # print(f"rco2_3", {rco2_3})
    # r_co2_part1 = rco2_1    
    rco2 = kn * (rco2_1 * rco2_2 - rco2_3)
    # print(f"rco2", {rco2})
    # r_co2_part1 =(Rg * T * c_co2 * ((1 - ((q / q_s0) ** (t_var))) ** (1 / t_var)) - q / (b_var * q_s0))
    # print(f"rco2_part1",{r_co2_part1})
    # r_co2 = kn *r_co2_part1
    # print(f'r_co2, {r_co2}')
    # print(f"b_var * qs_var, {rco2_term5}.")
    return rco2

"""
    Defines the differential equations for odes of DAC
    Arguments:
        y :  vector of the state variables:
                  y = all 15 states
        t :  time
        params :  vector of the parameters:
                  params = [V, r, T, c_co2_0, episl_r, v0]
    
"""

def deriv1(t, y, params):
    T_n, co2_n, q_n, T_n2, co2_n2, q_n2,T_n3, co2_n3, q_n3, T_n4, co2_n4, q_n4,T_n5, co2_n5, q_n5 = y # the rest of 12 vars a_n are not used, only for the success of solve_ivp
    V, T, c_co2_0, episl_r, volumetric_flow = params

    ###############   -----  Parameters depend on input  -----  ###############
    r = cube(V/(20*math.pi))
    v0 = volumetric_flow / (math.pi *r*r )
    L = V / (math.pi * (r ** 2))
    deltZ = L / 5.0  # 5 boxes in total
    a_s = 2 / r
    theta = (1 - episl_r) * ps * Cp_s + episl_r * pg * Cp_g
    t_var = t_h(T)

    # print(f"t_var in deriv", {t_var}) compare with the one in rco2
    b_var = b(T)
   # T_n, co2_n, q_n, T_n2, co2_n2, q_n2, T_n3, co2_n3, q_n3, T_n4, co2_n4, q_n4, T_n5, co2_n5, q_n5 == y
    # rco2_ first, rate of generation
    T1 = ((-v0  * pg* Cp_g) / (theta * deltZ)) * T_n + v0  * pg* Cp_g* T0 / (theta * deltZ)  + (1 - episl_r) * ps * deltH_co2 * (
        R_co2(T_n, co2_n, q_n,  b_var, t_var))/theta + a_s * h * (Tw - T_n)/theta
    # print(f"T1", {T1})
    co2_1dot = -v0 / (episl_r * deltZ) * co2_n + (v0 / (episl_r * deltZ))* c_co2_0 - (
        R_co2(T_n, co2_n, q_n,  b_var, t_var)) * (1 - episl_r) * ps/episl_r
    q1dot = R_co2(T_n, co2_n, q_n,  b_var, t_var)
    # print(f"energy balance in T1", {ener_balan(v0, theta, deltZ)})
    
    T2 = (-v0  * pg* Cp_g / (theta * deltZ)) * T_n2 + (v0  * pg* Cp_g / (theta * deltZ)) * T1 + (1 - episl_r) * ps * deltH_co2 * (
        R_co2(T_n2, co2_n2, q_n2,  b_var, t_var))/theta + a_s * h * (Tw - T_n2)/theta
    # print(f"T2", {T2})
    co2_2dot = -v0 / (episl_r * deltZ) * co2_n2 + (v0 / (episl_r * deltZ) )* co2_1dot - (
        R_co2(T_n2, co2_n2, q_n2,  b_var, t_var)) * (1 - episl_r) * ps/episl_r
    q2dot = R_co2(T_n2, co2_n2, q_n2,  b_var, t_var)
    # print(f"energy balance in T1", {ener_balan(v0, theta, deltZ)})

    T3 = (-v0  * pg* Cp_g / (theta * deltZ) )* T_n3 + (v0  * pg* Cp_g / (theta * deltZ)) * T2 + (1 - episl_r) * ps * deltH_co2 * (
        R_co2(T_n3, co2_n3, q_n3,  b_var, t_var))/theta + a_s * h * (Tw - T_n3)/theta
    co2_3dot = -v0 / (episl_r * deltZ) * co2_n3 + v0 / (episl_r * deltZ)* co2_2dot - (
        R_co2(T_n3, co2_n3, q_n3,  b_var, t_var)) * (1 - episl_r) * ps/episl_r
    q3dot = R_co2(T_n3, co2_n3, q_n3,  b_var, t_var)

    T4 = (-v0  * pg* Cp_g / (theta * deltZ) )* T_n4 + (v0  * pg* Cp_g / (theta * deltZ)) * T3 + (1 - episl_r) * ps * deltH_co2 * (
        R_co2(T_n4, co2_n4, q_n4,  b_var, t_var))/theta + a_s * h * (Tw -  T_n4)/theta
    co2_4dot = -v0 / (episl_r * deltZ)* co2_n4 + v0 / (episl_r * deltZ)* co2_3dot - (
        R_co2(T_n4, co2_n4, q_n4,  b_var, t_var)) * (1 - episl_r) * ps/episl_r
    q4dot = R_co2(T_n4, co2_n4, q_n4,  b_var, t_var)

    T5 = (-v0  * pg* Cp_g / (theta * deltZ)) * T_n5 + (v0  * pg* Cp_g / (theta * deltZ))* T4 + (1 - episl_r) * ps * deltH_co2 * (
        R_co2(T_n5, co2_n5, q_n5,  b_var, t_var))/theta + a_s * h * (Tw - T_n5)/theta
    co2_5dot = -v0 / (episl_r * deltZ) * co2_n5 + v0 / (episl_r * deltZ)* co2_4dot - (
        R_co2(T_n5, co2_n5, q_n5,  b_var, t_var)) * (1 - episl_r) * ps/episl_r
    q5dot = R_co2(T_n5, co2_n5, q_n5,  b_var, t_var)

    # result = np.array([T1, T2, T3, T4, T5, co2_1dot, co2_2dot, co2_3dot, co2_4dot, co2_5dot, q1dot, q2dot, q3dot, q4dot, q5dot]).reshape(-1, 1)

    return [T1, co2_1dot, q1dot, T2, co2_2dot, q2dot, T3, co2_3dot, q3dot, T4, co2_4dot, q4dot, T5, co2_5dot, q5dot]

# ------------------ User generated - Slider initial value -------------- #
V = 200.0  # volume
T = 293.0 # +273 ambrient temperature
c_co2_0 = .016349 # mol/m^3    
episl_r = 0.3  # void
volumetric_flow = 5 # m^3/s

# air humidity 
# no radius and length, nonly nr *reed V, july 6th

# ------------------ Initial Conditions to set up solve_ivp -------------- #
t0, tf = 0.0, 100.0 # 3hrs
co2_initial = 0
q_init_cond = 0
# init_cond = [20.000, 0.000, 0.000]
init_cond = [T, co2_initial, q_init_cond, T, co2_initial,q_init_cond, T, co2_initial, q_init_cond, T, co2_initial, q_init_cond, T,co2_initial, q_init_cond]
# ,20.000, 0.000, 0.000,20.000, 0.000, 0.000,20.000, 0.000, 0.000,20.000, 0.000, 0.000
params = [V, T, c_co2_0, episl_r, volumetric_flow]
tspan = np.linspace(t0, tf, 5)
soln = solve_ivp(deriv1, (t0, tf), init_cond, args=(params,), t_eval = tspan)  # init_cond = (T, c_co2_0, q0)
# soln = solve_ivp(deriv1, (t0, tf), init_cond, args=(params,), method = "Radau", rtol = 1e-5, atol = 1e-8)  # init_cond = (T, c_co2_0, q0)
# deriv1([t0, tf], )
print(soln)