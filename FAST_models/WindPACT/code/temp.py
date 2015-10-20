""" just dicking around
"""
import scipy.signal

#k   = 5.6e9
#Ir  = 34.6e3
#Ig  = 53.036
tau = 0.02
Kp = 5.14
Kd = 0.08-Kp*tau
print(Kd)
Kp = 3.0
Kd = -.0228
print(Kp,Kp*tau+Kd)

#b = [12947.976878612717, 831907.514450867]
#k   = 4.8e8
#Ir  = 2.998e4
#Ig  = 53.036
#tau = 0.02
#
#
#
#a0 = k/Ir
#b3 = Ig
#b2 = 1
#b1 = (1+Ig/Ir)*k
#b0 = 1*k/Ir
#
#
#Kp = b[1]/a0
#Kd = b[0]/a0-Kp*tau
#
#
#b = [a0*(Kp*tau+Kd), a0*Kp]
#a = [tau*b3, b3+tau*b2, b2*tau+b1, b1+tau*b0, b0]
#
#z,p,k = scipy.signal.tf2zpk(b,a)
#
#for i in range(4):
#    z = -np.cos(np.angle(p[i]))
#    wn = np.abs(p[i])
#    print('{:.2f}, {:.2f}'.format(z,wn))