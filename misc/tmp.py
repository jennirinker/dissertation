"""
"""

MinGridW = 55.
DX = 4.

n_el = np.round(MinGridW/DX)+np.round(MinGridW/DX)%2
n_x  = n_el + 1
GridW = DX * n_el

print(n_x,GridW)

MinGridW = 77.
DX = 5.

n_el = np.round(MinGridW/DX)+np.round(MinGridW/DX)%2
n_x  = n_el + 1
GridW = DX * n_el

print(n_x,GridW)

MinGridW = 109.
DX = 7.5

n_el = np.round(MinGridW/DX)+np.round(MinGridW/DX)%2
n_x  = n_el + 1
GridW = DX * n_el

print(n_x,GridW)

MinGridW = 141.
DX = 10.

n_el = np.round(MinGridW/DX)+np.round(MinGridW/DX)%2
n_x  = n_el + 1
GridW = DX * n_el

print(n_x,GridW)