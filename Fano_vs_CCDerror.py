from fanoCalculus import fano
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy import stats
import time
import sys
# ===================================================================
# Montecarlo Simulation Fano factor V deposited energy V ccd Error
# ===================================================================

# ===================================================================
# Global Parameters
# ===================================================================

E_eh = 3.6
fanoFactor = 0.1

# ===================================================================
# Functions
# ===================================================================


def gaus(x,a,x0,sigma):
    print sigma
    return a*np.exp(-(x-x0)**2/(2*sigma**2))

def gausFit(bins,n,measure,mu):
    # mask=np.isnan()
    poptS,pcovS = curve_fit(gaus,bins[:-1],n,p0=[1,mu,(mu*0.1)**0.5])
    print np.std(n)
    plt.plot(bins[:-1],n)
    plt.plot(bins,gaus(bins,*poptS))
    plt.show()
    print poptS
    print
    print pcovS
    # exit()
    # print fano(np)
    # print fano(poptS[1],pcovS[1,1],poptS[2],pcovS[2,2])
def inputs():
    Egamma    = input("Photons Energy[keV]: ")*10**3
    print Egamma
    sigmaCCD   = input("CCD electronic Noise [# electrons]: ")
    return float(Egamma/E_eh), float(sigmaCCD)

    # print samples[i],measure[i]
    # time.sleep(1)
# plt.hist(measure,3500,normed=True)
# m, s = stats.norm.fit(measure)
# print m,s
# x = np.linspace(measure.min(),measure.max(),len(measure))
# pdf_g = stats.norm.pdf(x, m, s)
# plt.plot(x, pdf_g)
# plt.show()

def main():
    if len(sys.argv) != 3:
        mu, sigmaCCD = inputs()
        # inputs(mu,sigma)
    else:
        mu = float(sys.argv[1])*10**3/E_eh
        sigmaCCD = float(sys.argv[2])
    print mu
    #
    sigma = np.sqrt(mu*fanoFactor) # Pairs electron-hole distribution sigma
    samples = np.round(np.random.normal(mu,sigma,100000)) # Rounded normal distribution (Electrons integers)
    measure = []
    #
    for i in range(len(samples)):
        measure = np.append(measure,np.random.normal(samples[i],sigmaCCD,1))
    n, bins, patches = plt.hist(measure,bins=np.linspace(int(measure.min()),int(measure.max()),int(measure.max())-int(measure.min())+1),align='left	')
    #gausFit(bins,n,measure,mu)
    #n, min_max, mean, var, skew, kurt = stats.describe(measure)
    #std=np.sqrt(var)
    # R = stats.norm.interval(0.95,loc=mean,scale=std)
    # print R
    # print fano(np.mean(measure),0,np.std(measure),0)
    plt.show()
main()
