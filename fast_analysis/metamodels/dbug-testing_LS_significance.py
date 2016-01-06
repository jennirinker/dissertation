"""
Testing least squares regression for RSM fitting
"""
import sys
libpath = 'C:\\Users\\jrinker\\Documents\\GitHub\\dissertation'
if (libpath not in sys.path): sys.path.append(libpath)
    
import JR_Library.main as jr
import numpy as np
import statsmodels.api as sm
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# create fake data
X1, X2 = np.meshgrid(np.arange(4),np.arange(5))
Y = 0.2*X1 + 1.1*(X1**2) + 0.7*X1*X2 + 0.3*(X2**2) + 3*np.random.rand(*X1.shape)
Xinp = np.hstack((X1.reshape(X1.size,1),X2.reshape(X2.size,1)))

def reg_m(y, x):
    ones = np.ones(len(x[0]))
    X = sm.add_constant(np.column_stack((ones,x[0])))
#    X = sm.add_constant(np.column_stack((x[0], ones)))
    for ele in x[1:]:
        X = sm.add_constant(np.column_stack((X,ele)))
    results = sm.OLS(y, X).fit()
    return results

# try to fit data
p_i = [2,2]
y = Y.reshape(Y.size)
ps = jr.GetAllPowers(p_i)
Xv = jr.myvander(Xinp,ps)
#results = sm.OLS(y,Xv).fit()
#print(results.summary())
results = reg_m(y, Xv[:,1:].T)
#print(results.summary())

coeffs, ps = jr.polyregression(Xinp,y,p_i)
py_coeffs = results.params
pvals     = results.pvalues
alpha = 0.05
print('Significance for alpha = {:.3f}'.format(alpha))
print('-------------------------------')
print(' x1 power x2 power mycoeff pycoeff pvalue    sig?')
p_red = np.empty(ps.shape)
coeffs_red = np.empty(py_coeffs.shape)
n_sig = 0
for i in range(len(ps)):
    if pvals[i] < alpha:
        sig = 'SIGNIFICANT'
        p_red[n_sig,:] = ps[i]
        coeffs_red[n_sig] = py_coeffs[i]
        n_sig += 1
    else:
        sig = 'not sig'
    print('{:6.0f} {:8.0f} {:8.2f} {:7.2f} {:6.2f}    {:7s}'.format(ps[i,0],ps[i,1],\
            coeffs[i],py_coeffs[i],pvals[i],sig))
p_red = p_red[:n_sig,:]
coeffs_red = coeffs_red[:n_sig]

Y_red = np.zeros(X1.shape)
for i_p in range(p_red.shape[0]):
    Y_red += coeffs_red[i_p] * (X1 ** p_red[i_p,0] + X2 ** p_red[i_p,1])


# plot data and fit
yhat = np.dot(Xv,coeffs)
fig = plt.figure(1,figsize=(10,10))
plt.clf()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(X1,X2,Y)
ax.plot_wireframe(X1,X2,yhat.reshape(X1.shape),color='b',label='all coeffs')
ax.plot_wireframe(X1,X2,Y_red,color='r',label='just significant')
plt.legend(fontsize='large')
#
## check significance of values
#e = y - yhat
#vary = 1./(y.size - coeffs.size)*(np.dot(e.T,e))
#C = vary * np.linalg.inv(np.dot(X.T,X))
#for i in range(coeffs.size):
#    print(C[i,i])
#    t = coeffs[i]/np.sqrt(C[i,i])
##    print(ps[i],t)

