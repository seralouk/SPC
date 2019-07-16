from scipy.stats import norm
import numpy as np
from scipy import integrate

def computeK2(alpha, m, si,ui):
    m = float(m)
    oldk = 0.0
    step = 100.0
    
    while (step > 1e-04):
        k = oldk
        temp = 0.0
        while (temp <= alpha):
            k = k + step
            EV = [0.0]
            IFP = ComputeIFP(alpha,k,m,si,ui)
            
            for m1 in np.arange(m/20, m/5, m/10): # will take on 2 values: m/20, 3m/20
                for pi in np.arange(0.1,0.6,0.1):
                    delta = np.inf # max(scores)
                    ITP = ComputeITP(alpha,k,m,si,pi,delta,m1,ui)
                    #print(ITP)
                    EV.append(EVstrong2(m,m1,si,pi,IFP,ITP))
                    #print(EVstrong2(m,m1,si,pi,IFP,ITP))
            temp = np.max(EV)
        oldk = k - step
        step = step / 10.
    relax = k - step
    return relax



def ComputeIFP(alpha,k,m,si,ui):
    uj = norm.ppf(1.-(k*alpha/(m*si)))
    rho = 1./np.sqrt(si) # threshold for subset
    def f(x):
        return (1.-norm.cdf((ui-rho*x) / (np.sqrt(1.-rho**2))))*norm.pdf(x)
    IFP, quad_err = integrate.quad(f, uj, np.Inf)
    return IFP

def ComputeITP(alpha,k,m,si,pi,delta,m1,ui):
    uj = norm.ppf(1.-(k*alpha/(m*si)))
    rho = 1./np.sqrt(si) # threshold for subset
    def f2(x):
        return (1.- norm.cdf((ui-delta*pi/rho-rho*x) / (np.sqrt(1.-rho**2))))*norm.pdf(x)
    ITP, quad_err = integrate.quad(f2, uj, np.inf)
    return round(ITP,5)
 
def EVstrong2(m,m1,si,pi,IFP,ITP):
    m0 = m - m1
    ev = m0 * si * IFP + m1 * (1.- pi) * si * ITP
    return round(ev,5)
