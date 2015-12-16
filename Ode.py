# import scipy
# import Odes
import numpy as np
# import pyximport; pyximport.install(setup_args={'include_dirs': np.get_include()})
# from Odessc import Hdt
from Es import H
from Jac import jac
from Parameters import U00, J0, dU0, mu
from PyParameters import tau
from scipy.integrate import ode
from math import sqrt
from itertools import chain
from timeit import timeit

L = 5
nmax = 7
dim = nmax+1

def f(t, y, *f_args):
    return 1j * y[0]

# f0tmp = [[i*(1+1j)/(2*L*dim)]*dim for i in range(L)]
# f0 = list(chain.from_iterable(f0tmp))#[item for sublist in f0tmp for item in sublist]
f0 = np.array([(i/dim)*(1+1j)/(2*L*dim) for i in range(L*dim)])

E0 = H(f0)
print E0
quit()

r = ode(Hdt).set_integrator('zvode', method='bdf')
r.set_initial_value(f0).set_f_params(2.0).set_jac_params(2.0)
t1 = 2*tau

def integ():
    r.integrate(t1)
print timeit(integ, number=1)

# r.integrate(t1)
print r.successful()

ff = r.y
Ef = H(ff)
print Ef

Q = Ef - E0
print Q

print ff
