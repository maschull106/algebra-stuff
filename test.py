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

#J = R.ideal(x**2, x*y**2, x*y*z, x*z**2, y**2*z**2, y*z**3, z**4, y**3-x*z)

f=GroebnerPolynomial.make(x**3+5*x*y**2+2*y*z+z**4+z-9, symbols=R.symbols)

I = R.ideal(x**2, x*y**2, x*y*z, x*z**2, y**2*z**2, y*z**3, z**4, y**3-x*z)
# I = R.ideal(x**2, y**2, z**2)
# I = R.ideal(x, y, z**2)
S = Quotient(I**2)
A = RingQuotient(S, I)
B = IdealQuotient(S, I)
l = S.basis


s = "x2 x3 x3y x3z x2y x2z xy2 x2y2 x2y2z xy3 xy4 xy3z xy2z xyz x2yz x2yz2 xyz2 xz2 x2z2 x2z3 xz3 y2z2 xy2z2 xy2z3 y3z2 y4z2 y3z3 y2z3 y2z4 y2z5 yz3 xyz3 xyz4 yz4 yz5 yz6 z4 xz4 xz5 z5 z6 z7 y3-xz y4-xyz y5-xy2z y5z-xy2z2 y4z-xyz2 y3z-xz2"


def make_poly(s: str):
    symbols = {'x': x, 'y': y, 'z': z}
    operators = {'+': 1, '-': -1}
    expression = []
    symb = 1
    exp = 1
    p = GroebnerPolynomial.make(1, order=R.order, symbols=R.symbols)
    for c in s:
        if c.isnumeric():
            if symb == 1:
                p *= float(c)
            else:
                exp = int(c)
        elif c in symbols:
            p *= symb**exp
            symb = symbols[c]
            exp = 1
        elif c in operators:
            p *= symb**exp
            expression.append(p)
            p = GroebnerPolynomial.make(operators[c], order=R.order, symbols=R.symbols)
            symb = 1
            exp = 1
    p *= symb**exp
    return sum(expression, start=p)


s = list(map(make_poly, s.split()))
s = R.sort_list(s)
