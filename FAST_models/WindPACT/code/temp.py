""" just dicking around
"""
import scipy.signal

#print('\n1.5 MW CertTest')
#
#ktor = 5.6e9
#Irot = 2962443.750
#Igen = 53.036
#GBR = 87.965
#
#wn = np.sqrt(ktor/Irot + ktor/(Igen*GBR**2))
#
#print(wn)
#print(wn/2/np.pi)
#
#print('\n5.0 MW Reference')
#
#ktor = 867.637E6
#Irot = 38759236.000
#Igen = 534.116
#GBR = 97
#
#wn = np.sqrt(ktor/Irot + ktor/(Igen*GBR**2))
#
#print(wn)
#print(wn/2/np.pi)

print('\n0.75 MW Rinker')

ktor = 1.3e8
Irot = 665139
Igen = 16.651
GBR = 62.832

Ieff = (Irot *(Igen*GBR**2))/(Irot+Igen*GBR**2)

wn = np.sqrt(ktor/Ieff)

c = 2*0.05*wn*Ieff

print(wn)
print(wn/2/np.pi)
print(c)

#print('\n1.5 MW Rinker')
#
#ktor = 483129639.713
#Irot = 2953248.500
#Igen = 56.442
#GBR = 87.965
#
#Ieff = (Irot *(Igen*GBR**2))/(Irot+Igen*GBR**2)
#
#wn = np.sqrt(ktor/Ieff)
#
#c = 2*0.05*wn*Ieff
#
#print(wn)
#print(wn/2/np.pi)
#print(c)

#print('\n5.0 MW Rinker')
#
#ktor = 2.3E9
#Irot = 64807728.000
#Igen = 438.855
#GBR = 160.85
#
#wn = np.sqrt(ktor/Irot + ktor/(Igen*GBR**2))
#
#print(wn)
#print(wn/2/np.pi)