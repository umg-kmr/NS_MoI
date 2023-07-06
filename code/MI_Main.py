import time
import numpy as np
import params as pp
import MI_Driver as X
import multiprocessing
import os

edat = pp.eos_name+'/'
eosnm = pp.eos_name.lower()+'.dat'

#Initial Conditions
V0 = 2*pp.a0*pp.r_core
H0 = pp.a0*(pp.r_core**2)


t1 = time.time() #To obtain the execution time
rhos = np.arange(pp.min_cd,pp.max_cd,pp.steps) 

if not os.path.exists('../output/'+edat):
   os.makedirs('../output/'+edat)

fname = '../output/'+edat+'MoIL.dat' #Path where the outputs will be stored
f = open(fname,'w')

#Function to implement multiprocessing and use the driver module
def task(m):
    r,P,M,rho,V,H,J = X.rk4_driver(pp.h_start,m*1e17,V0,H0)
    #Write even core densities to file
    #if iter==7:
     #   X.write_PMR(r,M,P,rho,m)
      #  iter = 10
   # elif iter%10==0 and iter!=0:
    #    X.write_PMR(r,M,P,rho,m)
    #iter+=1
    C = X.yC(r,V,H,M)[1]
    k2v = X.k2(r,V,H,M)
    l = X.lmbd(r,V,H,M)
    L = X.LMBD(r,V,H,M)

    return m,M[-1]/pp.Ms,r[-1]/1e3,C,k2v,l,L,X.I(J[-1],r[-1])
   
if __name__ == "__main__":
   with multiprocessing.Pool() as pool:
      res = pool.map(task,rhos)
      for x in res:
         rho1 = float(x[0])
         M1 = float(x[1])
         r1 = float(x[2])
         L1 = float(x[6])
         I1 = float(x[7])
         # Writes Mass (Solar Mass), Radius(Km), \Lambda, Moment of Inertia (kg m^2) to the file for the specified range of core densities.
         f.write(str(M1)+','+str(r1)+','+str(L1)+','+str(I1)+'\n')
         
f.close()
t2 = time.time()
t = (t2-t1)/3600

print(t,'Hours')
