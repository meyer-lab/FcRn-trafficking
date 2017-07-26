import sympy as sp

t, Qu, releaseF, Vp, sortF, Ve = sp.symbols('t Qu releaseF Vp sortF Ve', real=True, positive=True)
Cc = sp.Function("Cc")
Ce = sp.Function("Ce")

Cc_ = sp.Eq(Cc(t).diff(t), (Ce(t) * sortF * releaseF - Cc(t)) * Qu / Vp)
Ce_ = sp.Eq(Ce(t).diff(t), Qu / Ve * (Cc(t) + Ce(t) * ((1 - releaseF) * sortF - 1)))

eq = (Cc_, Ce_)

print(sp.classify_ode(Ce_))

output = sp.dsolve(eq, ics={Ce(0): 0})


print(output)
