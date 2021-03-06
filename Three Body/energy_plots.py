import numpy as np
import pandas as pd
import matplotlib.pyplot as mplt


read_data_e = np.genfromtxt('Lines_euler.txt',delimiter=',')
read_data_r = np.genfromtxt('Lines_rk4.txt',delimiter=',')
read_data_l = np.genfromtxt('Lines_le.txt',delimiter=',')

##print(read_data_e)

M1 = read_data_e[0,0] ; R1 = read_data_e[0,3]
M2 = read_data_e[0,1] ; R2 = read_data_e[0,4]
M3 = read_data_e[0,2] ; R3 = read_data_e[0,5]

x_cm_e = read_data_e[0,6] ; y_cm = read_data_e[0,7]

## X,Y Positions

X1_e = read_data_e[1:,0] ; Y1_e = read_data_e[1:,1]
X2_e = read_data_e[1:,2] ; Y2_e = read_data_e[1:,3]
X3_e = read_data_e[1:,4] ; Y3_e = read_data_e[1:,5]


S3 = np.sqrt( (X1_e-X2_e)**2 + (Y1_e-Y2_e)**2 )
S2 = np.sqrt( (X1_e-X3_e)**2 + (Y1_e-Y3_e)**2 )
S1 = np.sqrt( (X3_e-X2_e)**2 + (Y3_e-Y2_e)**2 )

## X,Y Momenta

KX1 = read_data_e[1:,6] ; KY1 = read_data_e[1:,7]
KX2 = read_data_e[1:,8] ; KY2 = read_data_e[1:,9]
KX3 = read_data_e[1:,10] ; KY3 = read_data_e[1:,11]


t = read_data_e[1:,12]

# Energy plots

K1_e = 0.5/M1*KX1**2 + 0.5/M1*KY1**2
K2_e = 0.5/M2*KX2**2 + 0.5/M2*KY2**2
K3_e = 0.5/M3*KX3**2 + 0.5/M3*KY3**2

G = 39.43279791722677
U_e= -G* ( M1*M2/S3 + M1*M3/S2 + M2*M3/S1)


M1 = read_data_r[0,0] ; R1 = read_data_r[0,3]
M2 = read_data_r[0,1] ; R2 = read_data_r[0,4]
M3 = read_data_r[0,2] ; R3 = read_data_r[0,5]

x_cm = read_data_r[0,6] ; y_cm = read_data_r[0,7]

## X,Y Positions

X1 = read_data_r[1:,0] ; Y1 = read_data_r[1:,1]
X2 = read_data_r[1:,2] ; Y2 = read_data_r[1:,3]
X3 = read_data_r[1:,4] ; Y3 = read_data_r[1:,5]


S3 = np.sqrt( (X1-X2)**2 + (Y1-Y2)**2 )
S2 = np.sqrt( (X1-X3)**2 + (Y1-Y3)**2 )
S1 = np.sqrt( (X3-X2)**2 + (Y3-Y2)**2 )

## X,Y Momenta

KX1 = read_data_r[1:,6] ; KY1 = read_data_r[1:,7]
KX2 = read_data_r[1:,8] ; KY2 = read_data_r[1:,9]
KX3 = read_data_r[1:,10] ; KY3 = read_data_r[1:,11]

U_r= -G* ( M1*M2/S3 + M1*M3/S2 + M2*M3/S1)
# Time

t = read_data_l[1:,12]

# Energy plots

K1_r = 0.5/M1*KX1**2 + 0.5/M1*KY1**2
K2_r = 0.5/M2*KX2**2 + 0.5/M2*KY2**2
K3_r = 0.5/M3*KX3**2 + 0.5/M3*KY3**2


M1 = read_data_l[0,0] ; R1 = read_data_l[0,3]
M2 = read_data_l[0,1] ; R2 = read_data_l[0,4]
M3 = read_data_l[0,2] ; R3 = read_data_l[0,5]

x_cm = read_data_l[0,6] ; y_cm = read_data_l[0,7]

## X,Y Positions

X1 = read_data_l[1:,0] ; Y1 = read_data_l[1:,1]
X2 = read_data_l[1:,2] ; Y2 = read_data_l[1:,3]
X3 = read_data_l[1:,4] ; Y3 = read_data_l[1:,5]


S3 = np.sqrt( (X1-X2)**2 + (Y1-Y2)**2 )
S2 = np.sqrt( (X1-X3)**2 + (Y1-Y3)**2 )
S1 = np.sqrt( (X3-X2)**2 + (Y3-Y2)**2 )

## X,Y Momenta

KX1 = read_data_l[1:,6] ; KY1 = read_data_l[1:,7]
KX2 = read_data_l[1:,8] ; KY2 = read_data_l[1:,9]
KX3 = read_data_l[1:,10] ; KY3 = read_data_l[1:,11]
U_l= -G* ( M1*M2/S3 + M1*M3/S2 + M2*M3/S1)
# Time

t = read_data_l[1:,12]

# Energy plots

K1_l = 0.5/M1*KX1**2 + 0.5/M1*KY1**2
K2_l = 0.5/M2*KX2**2 + 0.5/M2*KY2**2
K3_l = 0.5/M3*KX3**2 + 0.5/M3*KY3**2



mplt.plot(t,K1_e, color = "k", label = "Kinetic energy euler")
mplt.plot(t,K2_e, color = "k")
mplt.plot(t,K3_e, color = "k")

mplt.plot(t,K1_r, color = "g", label = "Kinetic energy rk4")
mplt.plot(t,K2_r, color = "g")
mplt.plot(t,K3_r, color = "g")


mplt.plot(t,K1_l, color = "b",linestyle='--', label = "Kinetic energy leapfrog")
mplt.plot(t,K2_l, color = "b",linestyle='--')
mplt.plot(t,K3_l, color = "b",linestyle='--')
mplt.xlabel("Time")
mplt.ylabel("Energy")

mplt.legend()
mplt.show()

mplt.plot(t, U_e,color ='k', label = "Potential energy euler")
mplt.plot(t, U_r,color = 'g', label = "Potential energy rk4")
mplt.plot(t, U_l,color ='b',linestyle='--', label = "Potential energy leapfrog")
mplt.legend()

mplt.xlabel("Time")
mplt.ylabel("Energy")
mplt.show()
mplt.plot(t,K1_l+K2_l+K3_l+U_e, color = "k", label = "Total energy euler")
mplt.plot(t,K1_l+K2_l+K3_l+U_r, color = "g", label = "Total energy rk4")
mplt.plot(t,K1_l+K2_l+K3_l+U_l, color = "b",linestyle='--', label = "Total energy leapfrog")

mplt.xlabel("Time")
mplt.ylabel("Energy")

##mplt.set_title("Energy plots")

mplt.legend()

mplt.show()
