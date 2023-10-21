a=0.00001;
b=13.00001;
p=2.3564029;
tar = 1;
counter = 0;
y0 = [p;0;1;0];
while abs(tar)>1.e-4 && counter <= 1000
    [t,y]=ode45(@newtonshoot,a:0.2:b,y0);
    figure(1),plot(t,y(:,1))
    tar = abs(y(end,1));
    disp([counter p tar])
    figure(1),plot(t,y(:,1))
    p = p - (y(end,1)/y(end,3));
    y0 = [p;0;1;0];
    counter = counter+1;
end
function [dy] = newtonshoot(t,y)
    v=y(2);
    z=y(3);
    w=y(4);

    dy(1,1) = v;
    dy(2,1) = -(1/t)*v + y(1)-2*y(1)^3;
    dy(3,1)= w;
    dy(4,1) = (1-6*v^2)*z + (-1/t)*w;
end





