# Formula for Theorem 3.1 of Dusart "Explicit estimates of some functions over primes"
# Would be in cython but something is wrong pkg-config?
# solve that another day

# Off the shelf parameters

# Height of verified RH
# Gourdon 2004 (unpublished)
H = 2445999556030
# Platt
# H = 30610046000

# Height for which sum of imaginary part of zeros is explicitly computed
T_0 = 1132490.982
sum_bound = 11.6377324

# Zero free region of form Re(s) >= 1 - (R log(Im(s)))^(-1)
R = 5.573412

# Rosser bound on number of zeros up to height T
a_1 = 0.137
a_2 = 0.443
a_3 = 1.588

def r(T):
    return a_1*log(T) + a_2*log(log(T)) + a_3
def q(y):
    return (a_1*log(y) + a_2)/(y*log(y)*log(0.5*y/pi))

# Ramare bound on number of zeros up to height T with Re(s) >= sigma.
def  c(sigma):
    return log(1 + 4.9/H*(3*H)^(8*(1-sigma)/3)*log(H)^(5-2*sigma)) + 51.5*log(H)^2/H

# Modified Bessel function
def z(b, m):
    return 2 * sqrt(m*b/R)
def J(b, m):
    zed = z(b,m)
    # Y = H
    Lwr = 2*m/zed*log(H)
    var('t')
    f = exp(-0.5*zed*(t+1/t))
    S = numerical_integral(f, Lwr, +Infinity)[0]
    return 0.25 * zed / m * S

# Faber and Kadiri ingredients

def B_1(T):
    return sum_bound + (0.5/pi + q(T_0)) * (log(T/T_0)*log(sqrt(T*T_0)*0.5/pi)) + 2*r(T_0)/T_0
def B_2(m,T):
    return (0.5/pi + q(T))*((1+m*log(T*0.5/pi))/(m^2*T^m) - (1+m*log(H*0.5/pi))/(m^2*H^m)) + 2*r(T)/T^(m+1)
def B_3(m):
    return (0.5/pi + q(H))*(1+m*log(H*0.5/pi))/(m^2*H^m) + 2*r(H)/H^(m+1)
def B_4(m,sigma):
    return c(sigma)*(1+1/m)/H^m
def B_5(b_0, m, sigma):
    if b_0 < m*R*log(H)^2:
        return c(sigma)*(1 + 0.5*R/b_0*log(H)^2/(m*R/b_0*log(H)^2-1))*exp(-b_0/R/log(H))/H^m
    else:
        return c(sigma)*(exp(-b/R/log(H))/H^m + J(b_0, m))

def legendre(m):
    var('u')
    # Legendre polynomial
    P = 0
    for k in range(m + 1):
        B = binomial(m, k)
        P = P + B*B*(u+1)^k*(u-1)^(m-k)
    return P * 2^(-m)

def M(P, m, delta, a):
    # expecting P to be absolute value of the mth Legendre polynomial with u = 1-2*u.
    fact = 1
    for i in range(m+1, 2*m+2):
        fact = fact*i
    f = P*(delta*u+a)^(m+1)
    A = numerical_integral(f, 0, 1)[0]
    return (fact*A)

def eps(b_0, m, sigma, delta, T):
    # upper bound on |\Psi(x) - x|/x for x >= exp(b_0).
    if T >= H:
        print "T too large, substituting H."
        T = H
    # typos in coefm in Dusart. Compare to Faber and Kadiri Erratum (2018)
    coefm = 2/delta^m * (B_5(b_0, m, sigma) + B_3(m)*(exp((sigma-1)*b_0) + exp(-b_0*sigma)) + B_4(m,sigma)*exp(-b_0*(1 - 1/R/log(H))) + 0.5*B_2(m, T)*exp(-b_0/2))
    coef0 = B_1(T)*exp(-b_0/2) + 0.5*exp(-3*b_0)
    const = 0.5*delta + log(2*pi)*exp(-b_0)
    # Integrate for a = 1 and a = 1 - delta
    P = abs(legendre(m))(u = 1-2*u)
    Mm0 = M(P, m, delta, 1)
    M00 = M(1, 0, delta, 1)
    Mm1 = M(P, m, delta, 1-delta)
    M01 = M(1, 0, delta, 1-delta)
    # print Mm0, M00, Mm1, M01
    epsilon = max(Mm0*coefm + M00*coef0, Mm1*coefm + M01*coef0) + const
    return epsilon
