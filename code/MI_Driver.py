import numpy as np
import pandas as pd
from scipy.interpolate import CubicSpline
import params as pp

edat = pp.eos_name+'/'
eosnm = pp.eos_name.lower()+'.dat'

#Name of the EOS being used
file_name = '../EOS_Tables/'+eosnm

#Declaring the constants
G = pp.G
Ms = pp.Ms
c = pp.c
alpha = pp.alpha
a0 = pp.a0
r_core = pp.r_core
f = pp.f
mtol = pp.mtol
des_tol = pp.des_tol

#Importing data from EOS tables
def EOS_load():
    cols = ['Baryon_Density','Pressure','Energy_Density']
    df = pd.read_csv(file_name,sep='\t',names=cols) #Loading the EOS from file, check the file structure of the EOS file and change 'sep' accordingly

    #Converting to appropriate units and storing the data
    kgm = 1e3
    fctr = c**2

    P = np.zeros(0)
    rho = np.zeros(0)

    for i in df.Pressure:
        P = np.append(P,i*kgm*fctr)

    for j in df.Energy_Density:
        rho = np.append(rho,j*kgm)

    for i in range(1,len(P)-1):
        if P[i]<P[i-1]:
            P = np.delete(P,i)
            rho = np.delete(rho,i)

    return P,rho

#Converting data to log
def convert_log():
    P,rho = EOS_load()
    P_log = np.zeros(0)
    rho_log = np.zeros(0)

    for i in P:
        P_log = np.append(P_log,np.log10(i))

    for j in rho:
        rho_log = np.append(rho_log,np.log10(j))

    return P_log, rho_log

#Fitting Cubic Splines
P_log, rho_log = convert_log()
inv_eos = CubicSpline(P_log,rho_log)
eos = CubicSpline(rho_log,P_log)

#TOV Equations
def dpdr(P,M,rho,r):
    pbyr = -((G*M*rho)/(r**2)) * ( 1 + ( P/( rho*(c**2) ) ) ) * ( 1 + (4*np.pi*(r**3)*P)/(M*(c**2)) ) * ( (1 - (2*G*M)/(r*(c**2))) **(-1) )
    return pbyr

#Function that returns dm/dr
def dmdr(r,rho):

    mbyr = 4*np.pi*(r**2)*rho
    return mbyr

#Differential equations involved in calculating the tidal love number 
def dVdr(H0,V0,m,r,P,rho):
    derv = num_der(1e-10,rho,P)
    z = ( 1 - (alpha*(m/r)) )**(-1)
    x = (-V0 * ( (2/r) +  ( z*( (alpha*(m/(r**2))) + (2*np.pi*r*alpha * ((P/c**2) - rho) ) ) ) )) - (H0*( ((-6/(r**2))*z) + ( 2*np.pi*alpha*z*((5*rho) + ((9*P)/(c**2)) + ((rho+(P/(c**2)))/derv) ) ) - ( (z**2)* (-alpha*( (4*np.pi*r*rho) - (m/(r**2)) ))**(2))))
    return x

def dHdr(V0):
    return V0

#Moment of Inertia Calculation

#Function to implement Kerr Metric element
def L(M,r,J0):
    th = 0
    w = 2*np.pi*f
    J = J0*w
    a = J/(M*c)
    rs = (2*G*M) / (c**2)
    sig = (r**2) + ((a**2)*(np.cos(th)**2))
    delt = (r**2) + (a**2) - (r*rs)
    return sig/delt 
    
def dJdr(P,M,rho,r,J0):
    Jbyr = ((8*np.pi)/3) * (r**4) * (rho + (P/(c**2)) ) * L(M,r,J0)
    return Jbyr
 
#Calculates the moment of inertia
def I(J,R):
    return ( J / (1 + ( (2*G*J)/((R**3)*(c**2))  ) ) )

#Defining a function to calculate numerical derivative required in the calculation of dH/dr
def num_der(h,rho0,P0):
    d = (eos(np.log10(rho0)+h) - eos(np.log10(rho0)))/h
    dn = d*(P0/(rho0*(c**2)))
    return dn

#RK4 stepper without adaptive size control
def rk4_stepper_noadpt(h,rho0,M0,P0,r0,H0,V0,J0):
    kp1 = h*dpdr(P0,M0,rho0,r0)
    km1 = h*dmdr(r0,rho0)
    kV1 = h*dVdr(H0,V0,M0,r0,P0,rho0)
    kH1 = h*dHdr(V0)
    kj1 = h*dJdr(P0,M0,rho0,r0,J0)

    kp2 = h*dpdr((P0+(kp1/2)),(M0+(km1/2)),rho0,(r0+(h/2)))
    km2 = h*dmdr((r0+(h/2)),rho0)
    kV2 = h*dVdr( (H0+(kH1/2)) , (V0+(kV1/2)) ,(M0+(km1/2)),(r0+(h/2)),(P0+(kp1/2)),rho0)
    kH2 = h*dHdr(V0+(kV1/2))
    kj2 = h*dJdr((P0+(kp1/2)),(M0+(km1/2)),rho0,(r0+(h/2)),(J0+(kj1/2)))

    kp3 = h*dpdr((P0+(kp2/2)),(M0+(km2/2)),rho0,(r0+(h/2)))
    km3 = h*dmdr((r0+(h/2)),rho0)
    kV3 = h*dVdr( (H0+(kH2/2)) , (V0+(kV2/2)) ,(M0+(km2/2)),(r0+(h/2)),(P0+(kp2/2)),rho0)
    kH3 = h*dHdr(V0+(kV2/2))
    kj3 = h*dJdr((P0+(kp2/2)),(M0+(km2/2)),rho0,(r0+(h/2)),(J0+(kj2/2)))

    kp4 = h*dpdr((P0+kp3),(M0+km3),rho0,(r0+h))
    km4 = h*dmdr((r0+h),rho0)
    kV4 = h*dVdr( (H0+kH3) , (V0+kV3) ,(M0+km3),(r0+h),(P0+kp3),rho0)
    kH4 = h*dHdr(V0+kV3)
    kj4 = h*dJdr((P0+kp3),(M0+km3),rho0,(r0+h),(J0+kj3))


    P0 = P0 + ( (1/6)*(kp1 + (2*kp2) + (2*kp3) + kp4) )

    M0 = M0 + ( (1/6)*(km1 + (2*km2) + (2*km3) + km4) )

    rho0 = 10**(inv_eos(np.log10(P0)))

    V0 = V0 + ( (1/6)*(kV1 + (2*kV2) + (2*kV3) + kV4) )

    H0 = H0 + ( (1/6)*(kH1 + (2*kH2) + (2*kH3) + kH4) )
    
    J0 = J0 + ( (1/6)*(kj1 + (2*kj2) + (2*kj3) + kj4) )

    r0 = r0 + h

    return P0, M0, rho0, r0, V0, H0, J0

#Adaptive step-size control
def adapt(des_tol,h_old,rhoi,Vi,Hi,Mi,Pi,ri,Ji):
    h=h_old
    #Calculate half step

    h=h/2
    kp1 = h*dpdr(Pi,Mi,rhoi,ri)
    km1 = h*dmdr(ri,rhoi)
    kV1 = h*dVdr(Hi,Vi,Mi,ri,Pi,rhoi)
    kH1 = h*dHdr(Vi)
    kj1 = h*dJdr(Pi,Mi,rhoi,ri,Ji)

    kp2 = h*dpdr((Pi+(kp1/2)),(Mi+(km1/2)),rhoi,(ri+(h/2)))
    km2 = h*dmdr((ri+(h/2)),rhoi)
    kV2 = h*dVdr( (Hi+(kH1/2)) , (Vi+(kV1/2)) ,(Mi+(km1/2)),(ri+(h/2)),(Pi+(kp1/2)),rhoi)
    kH2 = h*dHdr(Vi+(kV1/2))
    kj2 = h*dJdr((Pi+(kp1/2)),(Mi+(km1/2)),rhoi,(ri+(h/2)),(Ji+(kj1/2)))

    kp3 = h*dpdr((Pi+(kp2/2)),(Mi+(km2/2)),rhoi,(ri+(h/2)))
    km3 = h*dmdr((ri+(h/2)),rhoi)
    kV3 = h*dVdr( (Hi+(kH2/2)) , (Vi+(kV2/2)) ,(Mi+(km2/2)),(ri+(h/2)),(Pi+(kp2/2)),rhoi)
    kH3 = h*dHdr(Vi+(kV2/2))
    kj3 = h*dJdr((Pi+(kp2/2)),(Mi+(km2/2)),rhoi,(ri+(h/2)),(Ji+(kj2/2)))

    kp4 = h*dpdr((Pi+kp3),(Mi+km3),rhoi,(ri+h))
    km4 = h*dmdr((ri+h),rhoi)
    kV4 = h*dVdr( (Hi+kH3) , (Vi+kV3) ,(Mi+km3),(ri+h),(Pi+kp3),rhoi)
    kH4 = h*dHdr(Vi+kV3)
    kj4 = h*dJdr((Pi+kp3),(Mi+km3),rhoi,(ri+h),(Ji+kj3))

    P1 = Pi + ( (1/6)*(kp1 + (2*kp2) + (2*kp3) + kp4) )
    M1 = Mi + ( (1/6)*(km1 + (2*km2) + (2*km3) + km4) )
    V1 = Vi + ( (1/6)*(kV1 + (2*kV2) + (2*kV3) + kV4) )
    H1 = Hi + ( (1/6)*(kH1 + (2*kH2) + (2*kH3) + kH4) )
    J1 = Ji + ( (1/6)*(kj1 + (2*kj2) + (2*kj3) + kj4) )

    rho1 = 10**(inv_eos(np.log10(Pi)))
    r1 = ri + h

    #Second half-step
    kp1 = h*dpdr(P1,M1,rho1,r1)
    km1 = h*dmdr(r1,rho1)
    kV1 = h*dVdr(H1,V1,M1,r1,P1,rho1)
    kH1 = h*dHdr(V1)
    kj1 = h*dJdr(P1,M1,rho1,r1,J1)

    kp2 = h*dpdr((P1+(kp1/2)),(M1+(km1/2)),rho1,(r1+(h/2)))
    km2 = h*dmdr((r1+(h/2)),rho1)
    kV2 = h*dVdr( (H1+(kH1/2)) , (V1+(kV1/2)) ,(M1+(km1/2)),(r1+(h/2)),(P1+(kp1/2)),rho1)
    kH2 = h*dHdr(V1+(kV1/2))
    kj2 = h*dJdr((P1+(kp1/2)),(M1+(km1/2)),rho1,(r1+(h/2)),(J1+(kj1/2)))

    kp3 = h*dpdr((P1+(kp2/2)),(M1+(km2/2)),rho1,(r1+(h/2)))
    km3 = h*dmdr((r1+(h/2)),rho1)
    kV3 = h*dVdr( (H1+(kH2/2)) , (V1+(kV2/2)) ,(M1+(km2/2)),(r1+(h/2)),(P1+(kp2/2)),rho1)
    kH3 = h*dHdr(V1+(kV2/2))
    kj3 = h*dJdr((P1+(kp2/2)),(M1+(km2/2)),rho1,(r1+(h/2)),(J1+(kj2/2)))

    kp4 = h*dpdr((P1+kp3),(M1+km3),rho1,(r1+h))
    km4 = h*dmdr((r1+h),rho1)
    kV4 = h*dVdr( (H1+kH3) , (V1+kV3) ,(M1+km3),(r1+h),(P1+kp3),rho1)
    kH4 = h*dHdr(V1+kV3)
    kj4 = h*dJdr((P1+kp3),(M1+km3),rho1,(r1+h),(J1+kj3))

    P2 = P1 + ( (1/6)*(kp1 + (2*kp2) + (2*kp3) + kp4) )
    M2 = M1 + ( (1/6)*(km1 + (2*km2) + (2*km3) + km4) )
    V2 = V1 + ( (1/6)*(kV1 + (2*kV2) + (2*kV3) + kV4) )
    H2 = H1 + ( (1/6)*(kH1 + (2*kH2) + (2*kH3) + kH4) )
    J2 = J1 + ( (1/6)*(kj1 + (2*kj2) + (2*kj3) + kj4) )
    rho2 = 10**(inv_eos(np.log10(P1)))
    r2 = r1 + h

    P0, M0, rho0, rdel, V0, H0, J0 = rk4_stepper_noadpt(2*h,rhoi,Mi,Pi,ri,Hi,Vi,Ji)

    #Calculate tolerance and step-size
    tol_P = abs(P2-P0)
    tol_M = abs(M2-M0)
    tol_V = abs(V2-V0)
    tol_H = abs(H2-H0)
    tol_J = abs(J2-J0)

    #Prevent division by zero
    if tol_P == 0 or tol_M == 0 or tol_V == 0 or tol_H == 0 or tol_J == 0:
        tol_P = tol_M = tol_V = tol_H = tol_J = mtol

    h_newP = h_old*(des_tol/tol_P)**(1/5)
    h_newM = h_old*(des_tol/tol_M)**(1/5)
    h_newV = h_old*(des_tol/tol_V)**(1/5)
    h_newH = h_old*(des_tol/tol_H)**(1/5)
    h_newJ = h_old*(des_tol/tol_J)**(1/5)

    if h_newP > h_newM and h_newP > h_newV and h_newP > h_newH and h_newP > h_newJ:
        h_new = h_newP
    elif h_newM > h_newP and h_newM > h_newV and h_newM > h_newH and h_newM > h_newJ:
        h_new = h_newM
    elif h_newV > h_newP and h_newV > h_newM and h_newV > h_newH and h_newV > h_newJ:
        h_new = h_newV
    elif h_newJ > h_newP and h_newJ > h_newM and h_newJ > h_newV and h_newJ > h_newH:
        h_new = h_newJ
    else:
        h_new = h_newH

    #Limit on minimum step-size
    if h_new < mtol:
        h_new = h_old

    return h_new

def rk4_stepper(h,rho0,M0,P0,r0,H0,V0,J0):
    h = adapt(des_tol,h,rho0,V0,H0,M0,P0,r0,J0)

    kp1 = h*dpdr(P0,M0,rho0,r0)
    km1 = h*dmdr(r0,rho0)
    kV1 = h*dVdr(H0,V0,M0,r0,P0,rho0)
    kH1 = h*dHdr(V0)
    kj1 = h*dJdr(P0,M0,rho0,r0,J0)

    kp2 = h*dpdr((P0+(kp1/2)),(M0+(km1/2)),rho0,(r0+(h/2)))
    km2 = h*dmdr((r0+(h/2)),rho0)
    kV2 = h*dVdr( (H0+(kH1/2)) , (V0+(kV1/2)) ,(M0+(km1/2)),(r0+(h/2)),(P0+(kp1/2)),rho0)
    kH2 = h*dHdr(V0+(kV1/2))
    kj2 = h*dJdr((P0+(kp1/2)),(M0+(km1/2)),rho0,(r0+(h/2)),(J0+(kj1/2)))

    kp3 = h*dpdr((P0+(kp2/2)),(M0+(km2/2)),rho0,(r0+(h/2)))
    km3 = h*dmdr((r0+(h/2)),rho0)
    kV3 = h*dVdr( (H0+(kH2/2)) , (V0+(kV2/2)) ,(M0+(km2/2)),(r0+(h/2)),(P0+(kp2/2)),rho0)
    kH3 = h*dHdr(V0+(kV2/2))
    kj3 = h*dJdr((P0+(kp2/2)),(M0+(km2/2)),rho0,(r0+(h/2)),(J0+(kj2/2)))

    kp4 = h*dpdr((P0+kp3),(M0+km3),rho0,(r0+h))
    km4 = h*dmdr((r0+h),rho0)
    kV4 = h*dVdr( (H0+kH3) , (V0+kV3) ,(M0+km3),(r0+h),(P0+kp3),rho0)
    kH4 = h*dHdr(V0+kV3)
    kj4 = h*dJdr((P0+kp3),(M0+km3),rho0,(r0+h),(J0+kj3))


    P0 = P0 + ( (1/6)*(kp1 + (2*kp2) + (2*kp3) + kp4) )

    M0 = M0 + ( (1/6)*(km1 + (2*km2) + (2*km3) + km4) )

    rho0 = 10**(inv_eos(np.log10(P0)))

    V0 = V0 + ( (1/6)*(kV1 + (2*kV2) + (2*kV3) + kV4) )

    H0 = H0 + ( (1/6)*(kH1 + (2*kH2) + (2*kH3) + kH4) )
    
    J0 = J0 + ( (1/6)*(kj1 + (2*kj2) + (2*kj3) + kj4) )

    r0 = r0 + h

    return P0, M0, rho0, r0, V0, H0, J0

#Function for RK4 implementation
def rk4_driver(h,rho0,V0,H0):
    r0 = r_core
    M0=(4/3)*np.pi*(r0**3)*rho0
    P0 = 10**eos(np.log10(rho0))
    J0 = (2/5)*M0*(r0**2)

    P_rk = np.zeros(0)
    rho_rk = np.zeros(0)
    M_rk = np.zeros(0)
    r_rk = np.zeros(0)
    V_rk = np.zeros(0)
    H_rk = np.zeros(0)
    J_rk = np.zeros(0)

    while P0>=0:
        P_rk = np.append(P_rk,P0)
        M_rk = np.append(M_rk,M0)
        rho_rk = np.append(rho_rk,rho0)
        r_rk = np.append(r_rk,r0)
        V_rk = np.append(V_rk,V0)
        H_rk = np.append(H_rk,H0)
        J_rk = np.append(J_rk,J0)
        P0, M0, rho0, r0, V0, H0, J0 = rk4_stepper(h,rho0,M0,P0,r0,H0,V0,J0)

    return r_rk,P_rk,M_rk,rho_rk,V_rk,H_rk,J_rk

#Compactness and y calculation
def yC(r,V,H,M):
    y = r[-1]*(V[-1]/H[-1])
    C = ( G / (c**2) )*(M[-1]/r[-1])
    return y,C

#Tidal love number calculation
def k2(r,V,H,M):
    y,C = yC(r,V,H,M)
    tm1 = ( (8*(C**5)) / 5 ) * ((1 - (2*C))**2) * (2 + ((2*C)*(y-1)) - y )
    tm3 = (2*C)*(6 - (3*y) + ( (3*C)*((5*y)-8) ))
    tm4 = (4*(C**3)) * (13 - (11*y) + (C*((3*y) -2)) + ((2*(C**2))*(1+y)) )
    tm5 = 3*((1-(2*C))**2) * (2 - y + ((2*C)*(y-1)))*np.log(1 - (2*C))
    tm2 = (tm3 + tm4 + tm5)**(-1)
    k2 = tm1 * tm2
    return k2

#Tidal deformability calculation
def lmbd(r,V,H,M):
    l = (2/3)*k2(r,V,H,M)*(r[-1]**5)*(G**(-1))
    return l

#Dimensionless Tidal Deformability Parameter
def LMBD(r,V,H,M):
    L = (2/3)*k2(r,V,H,M)*(yC(r,V,H,M)[1]**(-5))
    return L

#Function to write out intermediate stellar structure out to the file. Optional.
def write_PMR(r,M,P,rho,rh):
    fname = '../output/'+edat+'integration_values/rho'+str(round(rh))+'.dat'
    fo = open(fname,'w')
    j=0
    for i in r:
        fo.write(str(i)+','+str(M[j])+','+str(P[j])+','+str(rho[j])+'\n')
        j+=1
    fo.close()
    return 0
