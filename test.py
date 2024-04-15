from algebra_stuff.polynomial_ring import *
import time


def timer(f, *args):
    t0 = time.time()
    res = f(*args)
    t1 = time.time()
    print(f"Executed in {t1-t0} seconds")
    return res


R = PolyRing(n=3)
x, y, z = R.symbols
I = R.ideal(x**2, y**2, y+z)

J = R.ideal(x**2, x*y**2, x*y*z, x*z**2, y**2*z**2, y*z**3, z**4, y**3-x*z)

I = R.ideal(x**2, y**2, z**2)
f=GroebnerPolynomial.make(x**3+5*x*y**2+2*y*z+z**4+z-9, symbols=I.symbols)
I.to_quotient_basis(f)

I = R.ideal(x**2, x*y**2, x*y*z, x*z**2, y**2*z**2, y*z**3, z**4, y**3-x*z)
S = Quotient(I**2)
A = RingQuotient(S, I)
B = IdealQuotient(S, I)
l = S.basis
