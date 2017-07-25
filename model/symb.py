import sympy as sp

t, Q, Qu, sortF, Vp, releaseF, Ve = sp.symbols('t Q Qu sortF Vp releaseF Ve')
Cc = sp.Function("Cc")(t)
Cp = sp.Function("Cp")(t)
Ce = sp.Function("Ce")(t)


Cc_ = sp.Eq(Cc.diff(t), Q * (Cp - Cc))
Cp_ = sp.Eq(Cp.diff(t), (Q * Cc - Q * Cp - Qu * Cp + Qu * Ce * sortF * releaseF) / Vp)
Ce_ = sp.Eq(Ce.diff(t), Qu / Ve * (Cp + Ce * ((1 - releaseF) * sortF - 1)))

sp.dsolve((Cc_, Cp_, Ce_), (Cc, Cp, Ce))
