import sympy as sp

# qq is ((1 - releaseF) * sortF - 1)
# yy is sortF * releaseF

t, Q, Qu, qq, Vp, yy, Ve = sp.symbols('t Q Qu qq Vp yy Ve', real=True, positive=True)
Cc = sp.Function("Cc")
Cp = sp.Function("Cp")
Ce = sp.Function("Ce")


Cc_ = sp.Eq(Cc(t).diff(t), Q * (Cp(t) - Cc(t)))
Cp_ = sp.Eq(Cp(t).diff(t), (Q * Cc(t) - Q * Cp(t) - Qu * Cp(t) + Qu * Ce(t) * yy) / Vp)
Ce_ = sp.Eq(Ce(t).diff(t), Qu / Ve * (Cp(t) + Ce(t) * qq))

eq = (Cc_, Cp_, Ce_)

print(sp.classify_ode(Ce_))

sp.dsolve(eq, hint='1st_power_series', ics={Cp(0): 0, Ce(0): 0})
