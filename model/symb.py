from sympy import *

t, Q, Qu, sortF, Vp = symbols('t Q Qu sortF Vp')
Cc = Function("Cc")(t)
Cp = Function("Cp")(t)
Ce = Function("Ce")(t)
releaseF = Symbol('releaseF')
Ve = Symbol('Ve')


Cc_ = Eq(Cc.diff(t, t), Q * (Cp - Cc))
Cp_ = Eq(Cp.diff(t, t), (Q * Cc - Q * Cp - Qu * Cp + Qu * Ce * sortF * releaseF) / Vp)
Ce_ = Eq(Ce.diff(t, t), Qu / Ve * (Cp + Ce * ((1 - releaseF) * sortF - 1)))

dsolve((Cc_, Cp_, Ce_), (Cc, Cp, Ce))
