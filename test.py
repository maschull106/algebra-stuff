from algebra_stuff.polynomial_ring import *


R = PolyRing(n=3)
x, y, z = R.symbols
I = R.ideal(x**2, y**2, y+z)
