"""
Testing least squares regression for RSM fitting
"""
import sys
libpath = 'C:\\Users\\jrinker\\Documents\\GitHub\\dissertation'
if (libpath not in sys.path): sys.path.append(libpath)
    
import JR_Library.main as jr
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from sklearn.cross_decomposition import PLSRegression

# create fake data
X1, X2 = np.meshgrid(np.arange(8),np.arange(7))
Y = 1.1*(X1**2) + 0.7*X1*X2 + 0.3*(X2**2) + 8*np.random.rand(*X1.shape)
Xinp = np.hstack((X1.reshape(X1.size,1),X2.reshape(X2.size,1)))


## try to fit data
#p_i = [2,2]
y = Y.reshape(Y.size)
#ps = jr.GetAllPowers(p_i)
#Xv = jr.myvander(Xinp,ps)
#results = jr.OLSfit(Xv, y)
##results = sm.OLS(y, Xv).fit()

pls2 = PLSRegression(n_components=2)
test = pls2.fit(Xinp, y)
yhat = pls2.predict(Xinp)

# plot data and fit
fig = plt.figure(2,figsize=(10,10))
plt.clf()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(X1,X2,Y)
ax.plot_wireframe(X1,X2,yhat.reshape(X1.shape),color='b',label='all coeffs')
#plt.legend(fontsize='large')
#
## check significance of values
#e = y - yhat
#vary = 1./(y.size - coeffs.size)*(np.dot(e.T,e))
#C = vary * np.linalg.inv(np.dot(X.T,X))
#for i in range(coeffs.size):
#    print(C[i,i])
#    t = coeffs[i]/np.sqrt(C[i,i])
##    print(ps[i],t)

