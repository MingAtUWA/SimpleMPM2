import math

# Modified Cam Clay
# Given Su, get szz (vertical stress)

def szz_from_Su_NC(phi, lmd, kpa, Su):
    phi *= math.pi / 180.0
    M = 6.0 * math.sin(phi) / (3.0 - math.sin(phi))
    K = 1.0 - math.sin(phi)
    alpha = (1.0 + 2.0 * K) / 3.0
    belta = (1.0 - K) * (1.0 - K) / (M * M * alpha) + alpha
    szz = math.exp((lmd * math.log(2.0*Su/M) + math.log(2.0) * (lmd-kpa) - (lmd-kpa) * math.log(belta) - kpa * math.log(alpha)) / lmd)
    print("szz = %f" % szz)

def Su_from_szz_NC(phi, lmd, kap, szz):
    phi *= math.pi / 180.0
    M = 6.0 * math.sin(phi) / (3.0 - math.sin(phi))
    K = 1.0 - math.sin(phi)
    alpha = (1.0 + 2.0 * K) / 3.0
    belta = (1.0 - K) * (1.0 - K) / (M * M * alpha) + alpha
    Su = M * 0.5 * math.exp((math.log(2.0)*(kpa-lmd) + (lmd-kpa)*math.log(belta) + kpa*math.log(alpha) + lmd*math.log(szz)) / lmd)
    print("Su = %f" % Su)

if __name__ == "__main__":
    # friction angle (degree)
    phi = 23.5
    # lambda
    lmd = 0.205
    # kappa
    kpa = 0.044
    
    # required undrained shear strain
    #Su = 57.2e3
    #szz_from_Su_NC(phi, lmd, kpa, Su)
    
    szz = 200000.0
    Su_from_szz_NC(phi, lmd, kpa, szz)
    