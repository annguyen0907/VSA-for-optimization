dim = 3*3;
pop_size = 3;
a = -30; b = 30;
x = a + (b-a) * rand(pop_size, dim);
m = (a + b)/2 * ones(1,dim);
s = 9   ;

f = gauss_distribution(x, m, s);
plot(x,f,'.')
grid on
title('Bell Curve')
xlabel('Randomly produced numbers')
ylabel('Gauss Distribution') 

function f = gauss_distribution(x, mu, s)
p1 = -.5 * ((x - mu)/s) .^ 2;
p2 = (s * sqrt(2*pi));
f = exp(p1) ./ p2; 
end