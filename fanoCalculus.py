import sys

# ======================================================================================
# The Obtjective of this script is to calculate the fano factor with the uncertainty
# ======================================================================================
#
#
def fano(mu,Dmu,sigma,Dsigma):
    F = (sigma**2) / mu
    print F
    dF = ((2*sigma*Dsigma/mu)**2 + (Dmu* (sigma/mu)**2)**2)**(0.5)
    return F, dF
#
def inputs():
    mu    = input("Fano's mean: ")
    Dmu   = input("Fano's mean uncertainty: ")
    sigma = input("Fano's Sigma: ")
    Dsigma = input("Fano's Sigma uncertainty: ")
    return mu, Dmu, sigma, Dsigma
#
def main():
    if len(sys.argv) != 5:
        mu, Dmu, sigma, Dsigma = inputs()
        # inputs(mu,sigma)
    else:
        mu = float(sys.argv[1])
        Dmu = float(sys.argv[2])
        sigma = float(sys.argv[3])
        Dsigma = float(sys.argv[4])

    F, dF = fano(mu,Dmu,sigma,Dsigma)

    print mu, Dmu, sigma, Dsigma
    print
    print F, dF
    
if __name__ == '__main__':
    main()
