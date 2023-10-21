a= 0 ;
b = 20;
N = 1000;
h = (b-a)/N;
j = 1:N+1;
r = (j-1)*h;
%yn = 1.6*sech(x);
%yn = ones(1,N+1);
%yn = zeros(1,N+1);
%vn = 2*exp(-r.^2) ;
vn = 2*exp(-r.^2) - 2*r.^2.*exp(-r);
tole =1;
n =0;
M = 300;

while (n<M)&& (tole>1.e-6) %%-Number of iterations less than M and tolerance more than e^{-6)
super =1+0.5*h./r(2:end-1); %%length (n-2) superdiagonal of A (coeff (j+1))
main = -2+h^2*(-1+6*(vn(2:N).^2)); %%Length (n-1)diagonal of A(coeff (j))
sub = 1-0.5*h./r(3:end); %%length (n-2)subdiagonal of A (coeff (j-1))
RHS = -vn(1:end-2)+2* vn(2:end-1)-vn(3:end) ...
-0.5*h./r(2:end-1).*(vn(3:end)-vn(1:end-2))+h^2*vn(2:end-1)...
-2*h^2*(vn(2:end-1).^3); %%The righthand side of the equation
A = diag(main,0)+ diag(super(1:end-1),1)+ diag(sub(1:end-1),-1); %%Matlab puts the tridiagonal matrix together
%%using the superdiagonal, diagonal and subdiagonal above.
%% Change entries in A to take into account boundary conditions
A(1,1)=-2/(3) - 2/(3*r(2)) + h^2*6*vn(1)^2 -h^2;
A(1,2)= 2/(3) + 2/(3*r(2));
z = A \ RHS'; %%(Aw=dd') N-1 vector(the correction in Newtons Method)
z1 = (4/3)*z(2)-(1/3)*z(3);
vn(2:end-1)=vn(2:end-1)+z'; %%Use only interior points
vn(1)=vn(1)+z1; %%The iterative step for y(1)
vn(end)= 0; %%boundary condition
figure(1),plot(r,vn)
n=n+1; %% Increase n (iteration) by 1
disp([n, norm(z)])
tole = norm(z);
pause
end