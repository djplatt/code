# Emit a C++ `euler_factors` map (as in examples/ec_5077.a1.cpp) for an
# elliptic curve over Q.  Run with:  sage euler_factors.sage
#
# Good prime  p :  L_p(T) = 1 - a_p T + p T^2   ->  {1, -a_p, p}
# Bad  prime  p :  L_p(T) = 1 - a_p T           ->  {1, -a_p}
#                  (a_p = +1 split mult., -1 nonsplit mult., 0 additive)

E = EllipticCurve([1, -1, 1, -45, 141])   
bound =  1000                  # include primes p < bound

N = E.conductor()

print("map<int64_t, vector<int64_t>> euler_factors = {")
for p in primes(bound):
    ap = E.ap(p)
    coeffs = [1, -ap] if N % p == 0 else [1, -ap, p]
    print("  {%d, {%s}}," % (p, ", ".join(str(c) for c in coeffs)))
print("};")
