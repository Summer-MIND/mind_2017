function f = ffitexp(x,xdata)

f = x(1)+(x(2)-x(1))*(1-exp(-xdata./x(3)));