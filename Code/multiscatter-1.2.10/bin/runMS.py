from numpy import *
import os
cmd='./multiscatter2 >out'
os.system(cmd)
  
datin=loadtxt('out')
Z=[]
lamb=0.00857
dZ=log10(pi**4/1e18/lamb**4/4*0.93)*10
for i in range(datin.shape[0]):
    Z.append(log10(datin[i,4]+1e-18)*10-dZ)
   
print Z 
