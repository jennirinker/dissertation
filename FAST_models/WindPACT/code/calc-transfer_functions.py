"""
Figure out natural frequencies, damping, of each section
"""
import numpy as np


# WP0.75 parameters
k = 1.3E8
c = 2.8E5
Ir = 665139.
Ig_HSS = 16.651
GBR = 62.832
Omega_rated = 28.648

# get drivetrain natural frequency and damping
Ig = Ig_HSS*GBR**2
num = [1]
den = [Ig,c*(1+Ig/Ir),k*(1+Ig/Ir)]

zeros = np.roots(num)
poles = np.roots(den)
wn    = np.abs(poles)[0]
zeta  = -np.cos(np.angle(poles))[0]

# SYSTEM 1: closed-loop controller
# default parameters: Ki = 2.22,Kp = 5.14, Kd = 0.0286
#c0 = 4
#a0,a1 = 5,0.08
#b0,b1 = 1.0,0.01
#Ki = c0
#Kp = a0
#Kd = a1 - Kp*b1

#Ki,Kp,Kd = 2.22, 5.14, 0.0286
#Ki,Kp,Kd = 2.22,5.14, 0.0286
#
#num_PID = [Kd,Kp,Ki]
#den_PID = [0,1,0]
#
#zeros = np.roots(num_PID)
#poles = np.roots(den_PID)
#wns   = np.abs(zeros)
#zetas = -np.cos(np.angle(zeros))
#
#print('\nPID controller:')
#print('---------------------')
#print('fns:   {:5.2f}   {:5.2f}'.format(*(wns/2/np.pi)))
#print('      ({:5.2f}) ({:5.2f})'.format(*(wns*60/(2*np.pi))))
#print('zetas: {:5.2f}   {:5.2f}'.format(*(zetas)))
#
##Ki,Kp,Kd = 2.22, 5.14, 0.0286
##Ki,Kp,Kd = 0,0,0
#Ki,Kp,Kd = 2.22,5, 20
#
Kd = 0.0286
KpKd = 160.72
KiKd = 30.6

Kp = KpKd * Kd
Ki = KiKd * Kd

num_PID = [Kd,Kp,Ki]
den_PID = [0,1,0]

zeros = np.roots(num_PID)
poles = np.roots(den_PID)
wns   = np.abs(zeros)
zetas = -np.cos(np.angle(zeros))
#
print(Kd,Ki,Kp)
print(zeros)
#print('\nClosed-loop system:')
#print('---------------------')
#print('fns:   {:5.2f}   '.format(*(wns/2/np.pi)))
#print('      ({:5.2f}) '.format(*(wns*60/(2*np.pi))))
#print('zetas: {:5.2f}   '.format(*(zetas)))



wns = 2*np.pi*np.logspace(-3,3,100)

C = (Kd*(1j*wns**2) + Kp*(1j*wns) + Ki)/(1j*wns)

plt.figure(10)
plt.clf()

plt.subplot(2,1,1)
plt.loglog(wns/2/np.pi,np.abs(C))
plt.axis('tight')

plt.subplot(2,1,2)
plt.semilogx(wns/2/np.pi,np.angle(C))

## blade pitch actuator
#num = [0,0,144.]
#den = [1,19.2,144.]
#
#zeros = np.roots(num)
#poles = np.roots(den)
#wns   = np.abs(poles)
#zetas = -np.cos(np.angle(poles))
#
#print('\nBlade pitch actuator:')
#print('---------------------')
#print('fns:   {:5.2f}   {:5.2f}'.format(*(wns/2/np.pi)))
#print('      ({:5.2f}) ({:5.2f})'.format(*(wns*60/(2*np.pi))))
#print('zetas: {:5.2f} {:5.2f}'.format(*(zetas)))
#
#
# drive train
#Ig = Ig_HSS*GBR**2
#num = [1]
#den = [Ig,c*(1+Ig/Ir),k*(1+Ig/Ir)]
#
#zeros = np.roots(num)
#poles = np.roots(den)
#wns   = np.abs(poles)
#zetas = -np.cos(np.angle(poles))
#
#
#print('\nDrive train:')
#print('------------')
#print('fns:   {:5.2f}   {:5.2f}'.format(*(wns/2/np.pi)))
#print('      ({:5.2f}) ({:5.2f})'.format(*(wns*60/(2*np.pi))))
#print('zetas: {:5.2f} {:5.2f}'.format(*(zetas)))
