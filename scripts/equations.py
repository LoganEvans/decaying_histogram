from pprint import pprint
from sympy import solve, Eq, S, nsolve

c_1 = S('c_1')
c_2 = S('c_2')
c_3 = S('c_3')
c_l = S('c_l')
c_h = S('c_h')
mu_1 = S('mu_1')
mu_2 = S('mu_2')
mu_3 = S('mu_3')
mu_l = S('mu_l')
mu_h = S('mu_h')

equations = [
        Eq(c_2, c_l + c_l),
        Eq(c_2, 2 * c_l),
        #Eq(c_2, 2 * c_h),
        Eq((mu_1 * c_1 + mu_2 * c_2) / (c_1 + c_2),
           (mu_1 * c_1 + mu_l * c_l) / (c_1 + c_l)),
        Eq((mu_2 * c_2 + mu_3 * c_3) / (c_2 + c_3),
           (mu_h * c_l + mu_3 * c_3) / (c_l + c_3))]

print solve(equations, [c_l, mu_l, mu_h])

