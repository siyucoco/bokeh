import math


###############    User generated - Slider initial value   ############### 
V= 100
r = 5
T = 25
c_co2_0 = 5
episl_r = 0.3 
v0 = 2
Tw = -5 # room temperature

###############  ---        Static Parameters          ---  ############### 
b0 = 93.0 * (10**(-5))
deltH_0 = 95.3 # calculate b
T0 = 353.15 # temeperature 
t0 = .37 # heterogeneity constant, in paper is denoted as t_h0
alpha = 0.33
chi = 0 
q_s0 = 3.40
R = 8.314
kT = 3.5*(10**3) #calculate rA
ρs = 880
deltH_co2 = 75 # calculate temeprature change 

### For Equation 4 : Enegergy Ballance 
ρg = 1.87 # ?
h = 13.8
Cp_g = 37.55 # J/molK
Cp_s = 1580 # J/molK
theta = (1-episl_r) * ρs * Cp_s +  episl_r * ρg * Cp_g

###############    -----   Input Parameters  -----  ############### 
L = V / (math.pi * r**2)
deltZ = L / 20
p_co2 = R * T * c_co2_0
a_s = deltZ / r



# Equations are calclulated in order 
def b(T):
    return( b0 * math.exp( (deltH_0/ (R * T0) ) * (T0/T - 1) ) )

def t_h(T):
    return ( t0 + alpha * (1 - T0 / T) )

def q_s(T):
    q_s0 * math.exp( chi * (1 - T / T0))

def q_eq(p_co2, T, q_s, b, t_h):
    return (q_s * b * p_co2 ) / ( ( 1 + ( (b * p_co2) **(t_h)) )** (1/t_h) )

# will be able to get rco_2
def R_co2(c_co2_0, q, T): 
    return ( kT * ( R * T * c_co2_0 * ( (1- ( (q / q_s)^t_h) )**(1/t_h) ) - q/ (b*q_s)) )

# ODE Part 
# Equation 1
# def co2_1(R_co2, episl_r, deltZ):
#     r_co2 = R_co2(c_co2_0, episl_r, deltZ)
#     co2_1 = ( ( -v0/ (episl_r * deltZ) ) * c_co2_0 + ( v0 / (episl_r * deltZ) ) - r_co2 * (1-episl_r) * ρs)
#     return co2_1

# Equation 2
# trying to find T1  
# dT/dt = ( ( -v0  * ρg* Cp_g / (theta * deltZ) ) * T1 + ( v0  * ρg* Cp_g / (theta * deltZ) ) * T0 - (1-episl_r) * ρs) * (deltH_co2) *r_co2 + a_s*h(Tw-T0)


# Equation 3
# dq/dt = r_co2


#  Graph co2, T1, and rco2
