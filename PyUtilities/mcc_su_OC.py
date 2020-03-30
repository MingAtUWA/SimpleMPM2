import math

# vertical stress
szz = 22000.0
# required undrained shear strain
Su = 50000.0
# friction angle (degree)
phi = 23.5
# lambda
lmd = 0.205
# kappa
kpa = 0.044

# main()
phi *= math.pi / 180.0
M = 6.0 * math.sin(phi) / (3.0 - math.sin(phi))
K = 1.0 - math.sin(phi)
alpha = (1.0 + 2.0 * K) / 3.0
belta = (1.0 - K) * (1.0 - K) / (M * M * alpha) + alpha
# minimum undrained shear strain
Su_min = M * 0.5 * math.exp(((kpa - lmd)*math.log(2.0) + (lmd - kpa)*math.log(belta*szz) + kpa*math.log(alpha*szz)) / lmd)
print("sv = %f, sh = %f\nmin Su is %f" % (szz, K * szz, Su_min))
if (Su < Su_min * 0.9999):
    print("Su too small\n")
    exit()

# initial p
p_ini = alpha * szz
pc = math.exp((lmd*math.log(2.0*Su/M) + (lmd - kpa)*math.log(2.0) - kpa*math.log(p_ini)) / (lmd - kpa))
print("p = %f, Su need pc %f " % (p_ini, pc))
if (p_ini > pc * 0.5):
    print("contractive\n")
elif (p_ini < pc * 0.5):
    print("dilative\n")
else:
    print("critical state\n")
