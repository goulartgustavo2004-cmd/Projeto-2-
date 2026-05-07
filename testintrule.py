from PyFET.IntRule1D import *


def f(xi):
    return 1+xi+xi**2+xi**3+xi**4
intrule = IntRule1D(100)
x,w=intrule.rule()
print(f"{x =}")
print(f"{w =}")
sum=0
for i in range(len(x)):
    sum+=w[i]*f(x[i])
print(f"{sum =}")