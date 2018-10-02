function y = legendreMP(lambda,P2,P4)
int1=0;
int2=pi;


fun = @(theta) 2*pi*exp(lambda(1)*(3*cos(theta).^2 - 1)/2 + lambda(2)*(35*cos(theta).^4 - 30*cos(theta).^2 + 3)/8).*sin(theta);



fun1 = @(theta) fun(theta).*(3*cos(theta).^2 - 1)/2;
fun2 = @(theta) fun(theta).*(35*cos(theta).^4 - 30*cos(theta).^2 + 3)/8;


y =  (integral(fun1,int1,int2)./integral(fun,int1,int2) - P2).^2 + (integral(fun2,int1,int2)./integral(fun,int1,int2) - P4).^2;


integral(fun,int1,int2);
end

