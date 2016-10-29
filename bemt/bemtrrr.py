## BEMT - Inflow Distribution and including the Tip Losses

import math

def prandtl_tip_loss(r,theta):
   global N_b, sol, cl_a
   err = 1
   lam_old = (sol*cl_a/16)*(math.sqrt(1 + 32*theta*r/(sol*cl_a)) - 1)
   while(err > 1e-5):
        f = 0.5*N_b*(1-r)/lam_old
        F = (2/math.pi)*math.acos(math.exp(-f))
        lam_new = (sol*cl_a/(16*F))*(math.sqrt(1 + 32*F*theta*r/(sol*cl_a)) - 1)
        err = abs(lam_new - lam_old)
        lam_old = lam_new
   return(lam_old)
   
   
# Constant Parameters
sol = 0.1 #Solidity
cl_a = 5.45 #Cl_a
N_b = 4   #Number of Blades

#Initialising Parameters                                                           

del_r = 0.010
r = [] #Non - Dimensional Radius
i = del_r
while i<1:
   r = r + [i]
   r[-1] = round(r[-1],3)
   i = i + del_r
N = len(r)

# Computation for varying theta
for theta in range(1,17,2):
   th = []     # Theta - pitch angle distribution
   del_ct = [] # Coefficient of thrust
   lam = [] # Inflow distribution
   
   for i in range(N):
       th = th + [theta*math.pi/180]
       
   for i in range(N):
       lam_t = prandtl_tip_loss(r[i],th[i])
       lam = lam + [lam_t]
       del_ct_t = N_b*0.5*(sol*cl_a)*(th[i]*r[i]**2 - lam_t*r[i])
       del_ct = del_ct + [del_ct_t]
 
   # Writing to output file
   f = "bemt_output_" + str(theta) + ".txt"
   out = open(f,'w')
   out.write('#r \t Inflow distribution\t\tThrust coefficient\n')
   for i in range(N):
      out.write(str(r[i])+ '\t' + str(lam[i])+ '\t' + str(del_ct[i])+ '\n')
   out.close()

print(r)
