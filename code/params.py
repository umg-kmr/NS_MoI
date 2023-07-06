#Change Equation of State name to the required one
eos_name = 'ALF2'
#Frequency of rotation of the neutron star in Hz
f = 300
h_start = 100 #Starting step-size
mtol = 1  #Minimum Step-Size
des_tol = 1e-8 #Desired Tolerance

#Range of core densities to be explored in S.I. Units/1e17
min_cd = 3.3
max_cd = 55
steps = 0.1
###########################################


#Constants
G = 6.6743e-11 #Gravitational Constant in S.I. Units
c = 2.997e8  #Speed of light in m/s
Ms = 1.9891e+30 #Mass of the Sun in kg
alpha = (2*G)/(c**2) #Redefining 2G/c^2 as \alpha
a0 = 0.1 #arbitrary value
r_core = 1e-10 #Staring value for the radius
