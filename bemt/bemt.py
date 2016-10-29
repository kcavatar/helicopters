## BEMT - inflow distribution

import math
# Constant parameters
sol = 0.1 #Solidity
cl_a = 5.45 #Cl_a

#Initialising parameters
lam = [] # Inflow dis
del_r = 0.050
r = [] #
i = del_r/2
th = []     # Theta - pitch angle distribution
del_ct = [] # Coefficient of thrust

# Computation
while i<1:
    r = r + [i]
    r[-1] = round(r[-1],3)
    th = th + [10*math.pi/180]
    i = i + del_r
    
N = len(r)
for i in range(N):
    th[i] = 10*(math.pi/180)/r[i]
for i in range(N):
    lam_t = (sol*cl_a/16)*(math.sqrt(1 + 32*th[i]*r[i]/(sol*cl_a)) - 1)
    lam = lam + [lam_t]
    del_ct_t = 0.5*(sol*cl_a)*(th[i]*r[i]**2 - lam_t*r[i])*del_r
    del_ct = del_ct + [del_ct_t]
C_T = sum(del_ct)

# Writing to output file

out = open("bemt_output.txt",'w')
out.write('#r \t Inflow distribution\n')
for i in range(N):
    out.write(str(r[i])+ '\t' + str(lam[i])+'\n')
out.close()

print(th)
