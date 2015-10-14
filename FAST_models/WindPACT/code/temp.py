""" just dicking around
"""

Cpmax = 0.45
TSR = 7.0
R  = 35.0
rho  = 1.225
GBR = 87.965
GenSpeed = 1800.


alpha = (rho*np.pi**3*R**5*Cpmax)/(1800*GBR**3*TSR**3)

GenSpeeds = np.linspace(0,2100,420)

