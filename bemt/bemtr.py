## BEMT - inflow distribution

import math
# Constant parameters
sol = 0.1 #Solidity
cl_a = 5.45 #Cl_a
N_b = 4   #Number of Blades

#Initialising parameters

del_r = 0.050
r = [] #
i = del_r/2
while i<1:
   r = r + [i]
   r[-1] = round(r[-1],3)
   i = i + del_r
N = len(r)

# Computation for varying theta
for theta in range(1,17,2):
   th = []     # Theta - pitch angle distribution
   del_ct = [] # Coefficient of thrust
   lam = [] # Inflow dis
   
   for i in range(N):
       th = th + [theta*math.pi/180]
       
   for i in range(N):
       lam_t = (sol*cl_a/16)*(math.sqrt(1 + 32*th[i]*r[i]/(sol*cl_a)) - 1)
       lam = lam + [lam_t]
       del_ct_t = N_b*0.5*(sol*cl_a)*(th[i]*r[i]**2 - lam_t*r[i])*del_r
       del_ct = del_ct + [del_ct_t]
   C_T = sum(del_ct)
 
   # Writing to output file
   f = "bemt_output_" + str(theta) + ".txt"
   out = open(f,'w')
   out.write('#Non-dimensional_radius \t Inflow distribution\t\tThrust coefficient\n')
   for i in range(N):
      out.write(str(r[i])+ '\t' + str(lam[i])+ '\t' + str(del_ct[i])+ '\n')
   out.close()

print(r)
