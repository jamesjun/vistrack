function x = wrap(x)

x = mod(x, 2*pi); 
x(x>pi) = x(x>pi) - 2*pi;