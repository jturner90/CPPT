#!/usr/bin/env python

import vegas
import math



def PS(x, a):
    p= 1. - (a**2 + x[0]**2 - x[2]**2)*(x[1]**2 + x[0]**2 - x[3]**2)/(4.*a*x[1]*x[0]**2)
    return p

def core(x,a,b):
    from numpy import cosh, sinh
    X = a + x[1] + x[2] + x[3]
    numer = X*sinh(X)*b
    denom = (X**2 + b**2)**2 * cosh(a) * cosh(x[1]) * sinh(x[2]) * sinh(x[3])
    return numer/denom

def integrand(x,a=1,b=0.1, eta3=1, eta4=1):
    # Subsitution rules and Jacobian
    x3 = eta3*abs(   a-x[0]) + (eta3*abs(a    + x[0]) - eta3*abs(a   -x[0]))*x[2]
    x4 = eta4*abs(x[1]-x[0]) + (eta4*abs(x[1] + x[0]) - eta4*abs(x[1]-x[0]))*x[3]
    J3 = eta3*abs(a    + x[0]) - eta3*abs(a   -x[0])
    J4 = eta4*abs(x[1] + x[0]) - eta4*abs(x[1]-x[0])
    x[2] = x3
    x[3] = x4
    # Return value
    r= 1./(2.*a) * PS(x,a) * core(x,a,b) * x[1] * J3 * J4
    return r



import sys
a=float(sys.argv[1])
b=float(sys.argv[2])
outfile = sys.argv[3]


ETA = [ [1., 1.], [1., -1.], [-1., 1.], [-1., -1.] ]

NEVAL=1e4

temp =[]
for eta in ETA:
    integ = vegas.Integrator([[-5., 5.], [-5., 5.], [0., 1.], [0., 1.]])
    # Training
    integ(         lambda x:integrand(x, a=a, b=b, eta3=eta[0], eta4=eta[1]), nitn=10, neval=NEVAL, alpha=0.1)#,  adapt_to_errors=True, alpha=0.1)
    # The real thing
    result = integ(lambda x:integrand(x, a=a, b=b, eta3=eta[0], eta4=eta[1]), nitn=100, neval=NEVAL, alpha=0.1)#, adapt_to_errors=True, alpha=0.1)
    temp.append(result.mean)


with open(outfile, "w") as f:
    f.write("%f\t%f\t%f\n"%(a,b,sum(temp)))

