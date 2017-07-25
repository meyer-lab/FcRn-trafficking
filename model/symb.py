import sympy as sp

# qq is ((1 - releaseF) * sortF - 1)
# yy is sortF * releaseF

t, Q, Qu, qq, Vp, yy, Ve = sp.symbols('t Q Qu qq Vp yy Ve', real=True, positive=True)
Cc = sp.Function("Cc")(t)
Cp = sp.Function("Cp")(t)
Ce = sp.Function("Ce")(t)


Cc_ = sp.Eq(Cc.diff(t), Q * (Cp - Cc))
Cp_ = sp.Eq(Cp.diff(t), (Q * Cc - Q * Cp - Qu * Cp + Qu * Ce * yy) / Vp)
Ce_ = sp.Eq(Ce.diff(t), Qu / Ve * (Cp + Ce * qq))

eq = (Cc_, Cp_, Ce_)

print(sp.classify_ode(Ce_))

sp.dsolve(eq, hint='1st_linear_Integral')
