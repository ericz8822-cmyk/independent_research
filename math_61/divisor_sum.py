import math
import cmath
from fractions import Fraction
from functools import lru_cache
import sympy as sp

p_table = [
    1, 1, 2, 3, 5, 7, 11, 15, 22, 30,
    42, 56, 77, 101, 135, 176, 231, 297, 385, 490,
    627, 792, 1002, 1255, 1575, 1958, 2436, 3010, 3718, 4565,
    5604, 6842, 8349, 10143, 12310, 14883, 17977, 21637, 26015, 31185,
    37338, 44583, 53174, 63261, 75175, 89134, 105558, 124754, 147273, 173525,
    204226, 239943, 281589, 329931, 386155, 451276, 526823, 614154, 715220, 831820,
    966467, 1121505, 1300156, 1505499, 1741630, 2012558, 2323520, 2679689, 3087735, 3554345,
    4087968, 4697205, 5392783, 6185689, 7089500, 8118264, 9289091, 10619863, 12132164, 13848650,
]

x, k = sp.symbols('x k', positive=True)
f = sp.sinh(sp.pi / k * sp.sqrt(sp.Rational(2, 3) * (x - sp.Rational(1, 24)))) / sp.sqrt(x - sp.Rational(1, 24))
df = sp.diff(f, x)

# ---------------------------------------------
# --------------------------------------------- VERIFIED
def round_up(number, decimals=0):
    multiplier = 10**decimals
    return math.ceil(number * multiplier) / multiplier
# ---------------------------------------------
# --------------------------------------------- VERIFIED
def divisor_sum(k): 
    total = 0
    for i in range(1, math.isqrt(k) + 1):
        if k % i == 0:
            if i < k // i:
                total += (i + k//(i))
            else:
                total += i
    return total
# ---------------------------------------------
# --------------------------------------------- VERIFIED
@lru_cache(maxsize=10000)
def p_int(n):
    if n <= 79:
        return p_table(n)
    total = 0
    for i in range(1, n + 1):
        total += divisor_sum(i) * p_int(n - i)
    return total // n
#remember to use p_int.cache_clear() after the excecution
# ---------------------------------------------
# --------------------------------------------- VERIFIED
def sawtooth(x):
    return 0 if x == int(x) else (x - math.floor(x) - 0.5)
# ---------------------------------------------
# --------------------------------------------- VERIFIED
def dedekind_sum(h, k):
    return sum(
        sawtooth(Fraction(r,k)) * 
        sawtooth(Fraction(h*r,k)) 
        for r in range(1, k))
# ---------------------------------------------
# --------------------------------------------- VERIFIED
def A_k(k,n):
    total = 0j
    for h in range(k):
        if math.gcd(h, k) == 1:
            s = dedekind_sum(h,k)
            exponent = math.pi * 1j * float((s - Fraction(2*n*h,k)))
            total += cmath.exp(exponent)
    return total.real
# ---------------------------------------------
# --------------------------------------------- VERIFIED
def p_int_est(n, t):
    total = 0
    for i in range(1, t+1):
        total += df.subs({x:n, k:i}) * math.sqrt(i) * A_k(i,n)
    total /= (math.pi * math.sqrt(2))
    return round_up(total, 20)