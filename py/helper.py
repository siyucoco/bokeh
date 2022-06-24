import math
from scipy.integrate import solve_ivp
from bokeh.core.property.instance import Instance
from bokeh.io import save
from bokeh.layouts import column
from bokeh.model import Model
from bokeh.models import CustomJS, Slider, Callback
from bokeh.plotting import ColumnDataSource, figure, show

# three plots , co2 as y and z as x


###############    User generated - Slider initial value   ############### 
V= 100.0 # volume
r = 5.0
T = 25.0
c_co2_0 = 5.0 # concentration
episl_r = 0.3 # void
v0 = 2.0 # initial vilocity
Tw = -5.0 # room temperature

###############  ---          Static Parameters        ---  ############### 
b0 = 93.0 * (10**(-5))
deltH_0 = 95.3 # calculate b
T0 = 353.15 # temeperature 
t0 = .37 # heterogeneity constant, in paper is denoted as t_h0
alpha = 0.33
chi = 0.0 
q_s0 = 3.40
R = 8.314
kT = 3.5*(10**3) #calculate rA
ρs = 880.0
deltH_co2 = 75.0 # calculate temeprature change 

# ------------------ For Equation 4 : Enegergy Ballance  --------------
ρg = 1.87 # ?
h = 13.8
Cp_g = 37.55 # J/molK
Cp_s = 1580.0 # J/molK

###############   -----  Parameters depend on input  -----  ############### 
L = V / (math.pi * r**2)
deltZ = L / 5.0 # 5 boxes in total
p_co2 = R * T * c_co2_0
a_s = deltZ / r
theta = (1-episl_r) * ρs * Cp_s +  episl_r * ρg * Cp_g

############# initial condition 
# init_cond = (T0, c_co2_0, co2_ad)


# Equations are calclulated in order 
def b(T):
    b = ( b0 ** ( (deltH_0/ (R * T0) ) * (T0/T - 1) ) )
    return b

def t_h(T):
    return ( t0 + alpha * (1 - T0 / T) )

def q_s(T):
    return ( q_s0 ** ( chi * (1 - T / T0)) )

def q(p_co2, T):
    b_var = b(T)
    t_var = t_h(T)
    qs_var = q_s(T)
    q= (qs_var * b_var * p_co2 ) / ( ( 1 + ( (b_var * p_co2) **(t_var)) )** (1/t_var) )
    return b_var, t_var, qs_var, q

# Calculate rco2_n (not ode), solve after have co2_n
def R_co2(co2_n, T): 
    b_var, t_h_var, q_s_var, q_var =  q(p_co2, T)
    r_co2 = ( kT * ( R * T * co2_n * ( (1- ( (q_var / q_s_var)**t_h_var) )**(1/t_h_var) ) - q_var/ (b_var*q_s_var)) )
    return r_co2
# ODE Part 
# Repetitive shortcut

# Equation 2
# trying to find T1  
# dT/dt = ( ( -v0  * ρg* Cp_g / (theta * deltZ) ) * T1 + ( v0  * ρg* Cp_g / (theta * deltZ) ) * T0 - (1-episl_r) * ρs) * (deltH_co2) *r_co2 + a_s*h(Tw-T0)
# equation Energy Balance shortcut

ener_balan_part1 = v0  * ρg* Cp_g 
def ener_balan(theta, deltZ): # replace v0  * ρg* Cp_g / (theta * deltZ) 
    return(ener_balan_part1/ (theta*deltZ) )
def ener_balan2(episl_r):
    return( (1-episl_r) * ρs * deltH_co2)
def ener_balan3(a_s, Tw, T0):
    return (a_s * h *(Tw-T0))
# def T_n():
#     return(-ener_balan(theta, deltZ) * Tn + ener_balan(theta, deltZ) * T0 + ener_balan2(episl_r)*rco2_n + ener_balan3(a_s, Tw, T0))

# Equation 1 Mass Balance : find co2_n

def mass_balan(episl_r, deltZ):
    return ( v0/ (episl_r * deltZ) )
def masss_balan2(episl_r, ρs):
    return( (1-episl_r ) * ρs )
# def co2_n(r_co2, episl_r, deltZ):
#     return(-mass_balan(episl_r, deltZ) * co2_1 + mass_balan(episl_r, deltZ) * c_co2_0 - r_co2 * masss_balan2(episl_r, ρs))
# -mass_balan(episl_r, deltZ) * co2_n + mass_balan(episl_r, deltZ) * c_co2_0 - r_co2 * masss_balan2(episl_r, ρs)

def deriv(t, y):
    T_n, co2_n, q_n = y
    
    T1 = -ener_balan(theta, deltZ) * T_n + ener_balan(theta, deltZ) * T0 + ener_balan2(episl_r)* r_co2_1 + ener_balan3(a_s, Tw, T0)
    co2_1 = -mass_balan(episl_r, deltZ) * co2_n + mass_balan(episl_r, deltZ) * c_co2_0 - r_co2_1 * masss_balan2(episl_r, ρs)
    r_co2_1 = R_co2(co2_1, T1)
    q_1 = r_co2_1

    T2 = -ener_balan(theta, deltZ) * T_n + ener_balan(theta, deltZ) * T0 + ener_balan2(episl_r)* r_co2_2 + ener_balan3(a_s, Tw, T0)
    co2_2 = -mass_balan(episl_r, deltZ) * co2_n + mass_balan(episl_r, deltZ) * c_co2_0 - r_co2_2 * masss_balan2(episl_r, ρs)
    r_co2_2 = R_co2(co2_2, T2)
    q_2 = r_co2_2

    T3 = -ener_balan(theta, deltZ) * T_n + ener_balan(theta, deltZ) * T0 + ener_balan2(episl_r)* r_co2_3 + ener_balan3(a_s, Tw, T0)
    co2_3 = -mass_balan(episl_r, deltZ) * co2_n + mass_balan(episl_r, deltZ) * c_co2_0 - r_co2_3 * masss_balan2(episl_r, ρs)
    r_co2_3 = R_co2(co2_3, T2)
    q_3 = r_co2_3

    T4 = -ener_balan(theta, deltZ) * T_n + ener_balan(theta, deltZ) * T0 + ener_balan2(episl_r)* r_co2_4 + ener_balan3(a_s, Tw, T0)
    co2_4 = -mass_balan(episl_r, deltZ) * co2_n + mass_balan(episl_r, deltZ) * c_co2_0 - r_co2_4 * masss_balan2(episl_r, ρs)
    r_co2_4 = R_co2(co2_4, T2)
    q_4 = r_co2_4

    T5 = -ener_balan(theta, deltZ) * T_n + ener_balan(theta, deltZ) * T0 + ener_balan2(episl_r)* r_co2_5 + ener_balan3(a_s, Tw, T0)
    co2_5 = -mass_balan(episl_r, deltZ) * co2_n + mass_balan(episl_r, deltZ) * c_co2_0 - r_co2_5 * masss_balan2(episl_r, ρs)
    r_co2_5 = R_co2(co2_5, T2)
    q_5 = r_co2_5
    return T1, co2_1, q_1, T2, co2_2, q_2, T3, co2_3, q_3, T4, co2_4, q_4, T5, co2_5, q_5

t0, tf = 0, 10
y0 = 1, 0, 0
soln = solve_ivp(deriv, (t0, tf), y0)
print(soln)
# Equation 3
# dq/dt = r_co2


#  Graph co2_n, T_n, and rco2_n



