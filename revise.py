import numpy as np
import math
import matplotlib.pyplot as plt

#constants
E_C0 = 0.9344
KC = 0.29
ES_0 = 0.65
AA = 7750000
BB = 2.9
S_H = 0
G_0_STAR = 1944
CHI = 0.94
P_A = 101000
H = 0.0533
M = 1.5
A_C = 4.97
V = 0.3
M2 = 3.16
LANDA = 0.0191
XI = 0.7
XALPHA = 1e25

print("E_C0:", E_C0)
print("KC:", KC)
print("ES_0:", ES_0)
print("AA:", AA)
print("BB:", BB)
print("S_H:", S_H)
print("G_0_STAR:", G_0_STAR)
print("CHI:", CHI)
print("P_A:", P_A)
print("H:", H)
print("M:", M)
print("A_C:", A_C)
print("V:", V)
print("M2:", M2)
print("LANDA:", LANDA)
print("XI:", XI)
print("ES_0", ES_0)
input()

# basic parameter
ND=6
# Axial strain
ep_a=[0.0]*ND
# Axial strain increment
dep_a=[0.0]*ND
# Stress increment
d_sig=[0.0]*ND

# Corrected stress
f_sig=[0.0]*ND
# Output Stress
sgm=[0.0]*ND
# Volumetric strain
ep_v=[0.0] 
# Stain tensor in total
ep_s=[0.0] * ND
# Stain tensor increment in total
dep_s=[0.0] * ND
# Plastic volumetric stain
ep_vp = 0.0 
# Plastic volumetric stain increment
dep_vp=0.0
# Plastic shear stain 
ep_ss=0.0
# Plastic shear stain increment
dep_ss=0.0


#Each Inter for the strain increment in total 
ddep_s=[0.0] * ND

#initial stress 
sgm_i=[0.0]*ND

#initial stress 
seff=[0.0] * ND

#corrected stress 
sgm_f=[0.0]*ND

DIMP=[0.0]*ND

# Elastic matrix
dem = np.zeros((ND, ND))

# Corrected elastic matrix
de = np.zeros((ND, ND))
rg = np.zeros((ND, ND))


print("TYPE OF TEST:")
print("DRAINED TRIAIAL BY STRAIN")
input()
# VE1---axial strain; NINC---increment
VE1=-0.15
NINC=100
print("AXIAL STRAIN:", VE1)
print("INCREMENT NUMBER:", NINC)
input()

dep_a[0] = VE1/NINC

DIMP[0] = dep_a[0] 

# print the initial stress state
sinit = [-1000000, -1000000, -1000000, 0, 0, 0]

for i in range(ND):
    sgm[i] = sinit[i]
    ep_s[i] = 0.0
    f_sig[i] = 0.0


for i in range(ND):
    seff[i]=-sinit[i]

# initial parameters
s = [0.0]*ND  

p_0=(seff[0] + seff[1] + seff[2]) / 3.0   

for i in range(ND):
    s[i] = 0.0
for i in range(3):
    s[i] = seff[i] - p_0
for i in range(3, ND):
    s[i] = seff[i]
        
q_0 = s[0] ** 2 + s[1] ** 2 + s[2] ** 2
q_0 += 2 * s[3] ** 2 + 2 * s[4] ** 2 + 2 * s[5] ** 2
q_0 = math.sqrt(3 * q_0 / 2)    
    
p_c = AA * S_H ** BB
# initial hardening parameters

e=ES_0-(KC*S_H)/(1+S_H)
# initial void ratio

eta= q_0/(p_0+p_c)
# initial stress ratio

e_c=E_C0-LANDA*(p_0/P_A)**XI
# initial void ratio of critical state line 

sp=e/e_c
# initial state parameter

G_v= math.sqrt(p_0 * P_A)*G_0_STAR * (1 - S_H) ** (-CHI) * ((2.97 - e) ** 2 / (1 + e)) 
# initial shear modulus

K_v=(2.0*(1+V)*G_v)/(3.0*(1-2*V))
# initial bulk modulus

E_yon=2*G_v*(1+V)
# initial young's modulus

d=M*sp**M2-eta
# initial dilantancy equation

f=q_0-eta*(p_0+p_c)
# initial yield function

d_epsi=[0.0] * ND
# increment strain tensor in total 

df_sig=[0.0] * ND
# increment stress in total 

ep_si=[0.0] * ND
# strain tensor in total 

q_values = []
#varying deviatoric stress

p_values = []
#varying effective stress

# Increment number (100)
for Increment in range(100):
    for i in range(ND):
        for j in range(ND):
            rg[i, j] = 0.0
            de[i, j] = 0.0                
   
    young = E_yon
    poiss = V
   
    def Elastic_matrix(young, poiss): 
        for i in range(ND):
            for j in range(ND):     
                dem[i, j] = 0.0
        coef= young/ (1.0 + poiss) / (1.0 - 2.0 * poiss) 
        coef1 = coef * (1.0 - poiss)
        coef2 = coef * poiss  
        for i in range(3):
            for j in range(3):
                if i == j:
                    dem[i, j] = coef1
                else:
                    dem[i, j] = coef2   
        for i in range(3, ND):
            dem[i, i] = coef * (1.0 - 2.0 * poiss) / 2.0      
    
        return dem
     
    dem = Elastic_matrix(young, poiss)

   # initial elastic matrix   
    for i in range(ND):
        for j in range(ND):
            de[i, j] = dem[i, j]
    for i in range(ND):
        for j in range(ND):
            rg[i, j] = de[i, j]      
    
    rg[0, 0]+=XALPHA
    
    for i in range(ND):
        f_sig[i]+=d_sig[i]
        ep_si[i]=0.0
      
    sgm_x = [0.0] * ND
    for i in range(ND):
        sgm_x[i] = sgm[i] - sinit[i]
        df_sig[i] = 0.0   
      
    for i in range(ND):
        df_sig[i] = f_sig[i] - sgm_x[i]
   
    df_sig[0]=XALPHA * dep_a[0]

    dep_s = np.linalg.inv(rg).dot(df_sig) 
    
    for i in range(ND):
        ep_si[i] += dep_s[i]
      
    for i in range(ND):
        sgm_i[i]=sgm[i]

    young = E_yon
    poiss = V
    dem = Elastic_matrix(young, poiss) 
    
    dep_sp=[0.0]*ND
    class BreakOutOfFunction(Exception):
      pass 
    def update_equations():
        global ep_vp, dep_s, sgm_i, dep_sp, ep_ss, dep_ss, ep_v, dep_vp, ep_si,ddep_s, f, s, sgm_f, sgm_a, p, q,e, e_c, G_v, K_v, E_yon, sp, d, eta, p_c,de, e, A_hm, d_f_sgm, d_g_sgm, x_landa, ddeps_e, ddeps_p, dem
        
        # Each increment divided into 5 iteration     
        for i in range(ND):
         sgm_f[i]=-sgm_i[i]
         ddep_s[i]=-dep_s[i]/5
         dep_sp[i] = 0.0
        
        sgm_a=[0.0]* ND
        ddeps_p=[0.0]* ND        
        
        try: 
        # Inter number   
          for INTER in range(5):     
          
            E_0=ES_0-(KC*S_H)/(1+S_H)
            P_C0=AA*S_H*BB 
            
            INTER += 1
            print(INTER)
          
            for i in range(ND):
                 sgm_a[i]=sgm_f[i]
                 ddeps_p[i]=0.0

            for i in range(ND):
                 temp = 0
                 for j in range(ND):
                   temp += dem[i, j]*ddep_s[j]
                 sgm_f[i] = sgm_a[i] + temp
            
            # plastic behavior 
            for i in range(ND):
               for j in range(ND):
                 de[i, j] = dem[i, j]         
            
            p= np.mean(sgm_a[:3])     
            
            for i in range(ND):
                 s[i] = 0.0
            for i in range(3):
                 s[i] = sgm_a[i] - p
            for i in range(3, ND):
                 s[i] = sgm_a[i]
    
            q = s[0] ** 2 + s[1] ** 2 + s[2] ** 2
            q = q + 2 * s[3] ** 2 + 2 * s[4] ** 2 + 2 * s[5] ** 2
            q = math.sqrt(3.0 * q / 2.0)
    
            f=q-M*(p+p_c)  
            
            print(f)
            input(f)
            
            if f < 1e-3:
                for i in range(ND):
                  dep_sp[i] += ddeps_p[i]
                for i in range(ND):
                  sgm_f[i]=-sgm_f[i]
                  dep_sp[i]=-dep_sp[i]
                raise BreakOutOfFunction           
                       
            dfdq=1.0
         
            dfdp=-eta
  
            dgdq=1.0
        
            dgdp= d
        
            d_f_sgm=ND*[0.0]
        
            d_g_sgm=ND*[0.0]
        
            if q < 1e-6:
                xx = 0
            else: 
                xx = 3/(2*q)
            for i in range(3):
             d_f_sgm[i] = dfdp*(1 / 3) + dfdq*xx*s[i]
             d_g_sgm[i] = dgdp*(1 / 3) + dgdp*xx*s[i]

            for i in range(3, ND):
             d_f_sgm[i] = dfdq*xx*s[i]
             d_g_sgm[i] = dgdq*xx*s[i]         
                             
            if eta <= 1e-6:
              A_hm = 0
            else:
              A_hm = H*G_v*(M/eta-sp)-A_C*eta*P_C0 
                   
            ff=0.0
            xx=0.0
            for i in range(ND):
              for j in range(ND):
                 xx+=de[i,j]*ddep_s[j]
                 ff+=d_f_sgm[i]*xx   
           
            ss=0.0
            yy=0.0
            for i in range(ND):
              for j in range(ND):
                yy+=de[i, j]*d_g_sgm[j]
                ss+=d_f_sgm[i]*yy
        
            ss += A_hm
            if ss <= 1e-6:
                 ss = 0
                 x_landa=0
            else:
             x_landa= ff/ss         
                
            if x_landa<= 0.0:
               x_landa = 0        
       
            ddeps_p=ND*[0.0]
            ddeps_e=ND*[0.0]
        
            for i in range(ND):
                ddeps_p[i]+=x_landa * d_g_sgm[i]
                ddeps_e[i]=ddep_s[i]-ddeps_p[i]
            for i in range(ND):
                temp = 0
                for j in range(ND):
                    temp += de[i, j]*ddeps_e[j]
                    sgm_f[i] = sgm_a[i] + temp
         

            p = np.mean(sgm_f[:3]) 
                
            for i in range(ND):
              s[i] = 0.0
            for i in range(3):
              s[i] = sgm_f[i] - p
            for i in range(3, ND):
              s[i] = sgm_f[i]
    
            q = s[0] ** 2 + s[1] ** 2 + s[2] ** 2
            q = q + 2 * s[3] ** 2 + 2 * s[4] ** 2 + 2 * s[5] ** 2
            q = math.sqrt(3.0 * q / 2.0)
              
            for i in range(ND):
              ep_si[i]+=ddep_s[i]
             #corrected strain
         
            ep_v =np.sum(ep_si[:3])/ 3.0
            #corrected volumetric strain
                  
            e=E_0-(1+E_0)*ep_v 
            #corrected void ratio
       
            dep_ss = x_landa* dgdq
            #corrected shear strain increment
         
            ep_ss += dep_ss
            #corrected plastic shear strain 
            
            dep_vp = x_landa* dgdp  
            ep_vp += dep_vp
            #corrected plastic volumetric strain 
         
            p_c=P_C0-A_C*P_C0*dep_ss  
            #corrected hardening parameters
          
            eta=q/(p + p_c)    
            #corrected stress ratio
         
            e_c=E_C0-LANDA*(p/P_A)**XI
            #corrected void ratio of critical state
         
            sp=e/e_c
            #corrected state parameters
            
            G_v= G_0_STAR * (1 - S_H) ** (-CHI) * (2.97 - e) ** 2 / (1 + e) * math.sqrt(p * P_A)
            #corrected shear modulus
         
            K_v=2*(1+V)*G_v/(3*(1-2*V))
            #corrected bulk modulus
         
            E_yon=2*G_v*(1+V)
            #corrected young's modulus
         
            d=M*sp**M2-eta
            #corrected dilatancy equation
               
            young=E_yon
            poiss= V
         
            dem = Elastic_matrix(young, poiss)
            #corrected elastic matrix
            for i in range(ND):
             for j in range(ND):
                de[i, j] = dem[i, j]         
        
            while True:
             p = np.mean(sgm_f[:3])
          
             for i in range(ND):
                s[i] = 0.0
             for i in range(3):
                s[i] = sgm_f[i] - p
             for i in range(3, ND):
                s[i] = sgm_f[i]
    
             q = s[0] ** 2 + s[1] ** 2 + s[2] ** 2
             q = q + 2 * s[3] ** 2 + 2 * s[4] ** 2 + 2 * s[5] ** 2
             q = math.sqrt(3.0 * q / 2.0)      
            
             f=q-M*(p+p_c) 

             if f< 1e-3:
                for i in range(ND):
                  dep_sp[i] += ddeps_p[i]
                for i in range(ND):
                    sgm_f[i] = -sgm_f[i]
                    dep_sp[i] = -dep_sp[i] 
                raise BreakOutOfFunction
         
             dfdq=1.0      
             dfdp=-eta
             dgdq=1.0
             dgdp=d
         
             if q < 1e-6:
                xx = 0
             else: 
                xx = 3/(2*q)
             for i in range(3):
                d_f_sgm[i] = dfdp*(1 / 3) + dfdq*xx*s[i]
                d_g_sgm[i] = dgdp*(1 / 3) + dgdp*xx*s[i]

             for i in range(3, ND):
                d_f_sgm[i] = dfdq*xx*s[i]
                d_g_sgm[i] = dgdq*xx*s[i]         
         
             if eta <= 1e-6:
                A_hm = 0
             else:
                A_hm = H*G_v*(M/eta-sp)-A_C*eta*P_C0 
            
             ff=0.0
             xx=0.0
             for i in range(ND):
                for j in range(ND):
                    xx+=de[i,j]*ddep_s[j]
                    ff+=d_f_sgm[i]*xx   
             ss=0.0
             yy=0.0
             for i in range(ND):
                for j in range(ND):
                    yy+=de[i, j]*d_g_sgm[j]
                    ss+=d_f_sgm[i]*yy
    
             ss += A_hm
             if ss <= 1e-6:
                ss = 0
                x_landa=0
             else:
                x_landa= ff/ss         
                
             if x_landa<= 0.0:
                  x_landa = 0        
         
             for i in range(ND):
                ddeps_p[i]+=x_landa * d_g_sgm[i]
                ddeps_e[i]=ddep_s[i]-ddeps_p[i]
             for i in range(ND):
                temp = 0
                for j in range(ND):
                    temp += de[i, j]*ddeps_e[j]
                    sgm_f[i] = sgm_a[i] + temp
                      
         
             p = np.sum(sgm_f[:3]) / 3.0 
         
             for i in range(ND):
                s[i] = 0.0
             for i in range(3):
                s[i] = sgm_f[i] - p
             for i in range(3, ND):
                s[i] = sgm_f[i]
    
             q = s[0] ** 2 + s[1] ** 2 + s[2] ** 2
             q = q + 2 * s[3] ** 2 + 2 * s[4] ** 2 + 2 * s[5] ** 2
             q = math.sqrt(3.0 * q / 2.0)
           
             for i in range(ND):
                ep_si[i]+=ddep_s[i]
             #corrected strain
         
             ep_v =np.sum(ep_si[:3])/ 3.0
             #corrected volumetric strain
                  
             e=E_0-(1+E_0)*ep_v 
             #corrected void ratio
         
             dep_ss = x_landa* dgdq
             #corrected shear strain increment
         
             ep_ss += dep_ss
             #corrected plastic shear strain 
            
             dep_vp += x_landa* dgdp  
             ep_vp += dep_vp
             #corrected plastic volumetric strain 
         
             p_c=P_C0-A_C*P_C0*dep_ss  
             #corrected hardening parameters
          
             eta=q/(p + p_c)    
             #corrected stress ratio
         
             e_c=E_C0-LANDA*(p/P_A)**XI
             #corrected void ratio of critical state
         
             sp=e/e_c
             #corrected state parameters
                 
             G_v= G_0_STAR * (1 - S_H) ** (-CHI) * (2.97 - e) ** 2 / (1 + e) * math.sqrt(p * P_A)
             #corrected shear modulus
         
             K_v=2*(1+V)*G_v/(3*(1-2*V))
             #corrected bulk modulus
         
             E_yon=2*G_v*(1+V)
             #corrected young's modulus
         
             d=M*sp**M2-eta
             #corrected dilatancy equation
               
             young=E_yon
             poiss= V
         
             dem = Elastic_matrix(young, poiss)
             #corrected elastic matrix
             for i in range(ND):
                for j in range(ND):
                      de[i, j] = dem[i, j]  
        except BreakOutOfFunction:
          pass
    update_equations()
         
    dep_se=[0.0]*ND
    
    for i in range(ND):
        dep_se[i] = dep_s[i]-dep_sp[i] 
    
    d_s=[0.0]*ND

    for i in range(ND):
        temp = 0
        for j in range(ND):
            temp +=de[i, j]*dep_se[j]
        d_s[i]+= temp

    for i in range(ND):
        sgm_f[i]=sgm_i[i]+d_s[i]

    for i in range(ND):
        sgm_x[i]=sgm_f[i]-sinit[i]
      
    for i in range(ND):
        df_sig[i]=f_sig[i]-sgm_x[i]     

    for i in range(ND):
        sgm[i]=sgm_f[i]

    for i in range(ND):
        ep_s[i]+=ep_si[i]
    
    p = -np.sum(sgm[:3]) / 3.0 
    #the effective stress after each increment
    
    q=-(sgm[0]-sgm[2])
    #the deviatoric stress after each increment
    
    qa = 3.0 * p_c + 3000000.0
    #the critical state of deviatoric stress
       
    q_values.append(q)
    p_values.append(p)   
    
    print(f"Increment  {Increment + 1}: p = {p}, q = {q}")
    
    #if q greater than qa, than skip this increment
    
    Increment+=1
    print(Increment)
    
plt.plot(range(len(q_values)), q_values, label='q')
plt.plot(range(len(p_values)), p_values, label='p')
plt.xlabel('Iteration')
plt.ylabel('Value')
plt.legend()
plt.show()