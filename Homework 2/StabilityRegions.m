close all
clc
% X,Y data points in complex plane
x = -6:.01:2;
y = -4:.01:4;
[X,Y] = meshgrid(x,y);

% forward Euler
f = @(z) 1+z;
Z = 1*(abs(f(X+1i*Y))<1);
subplot(3,2,1)
contourf(X,Y,Z,[1 1])
title('Forward Euler','Interpreter','latex','FontSize',14)
xlabel('$Re(z)$','Interpreter','latex','FontSize',14)
ylabel('$Im(z)$','Interpreter','latex','FontSize',14)
grid on,xline(0),yline(0),axis square
subplot(3,2,6)
contour(X,Y,Z,[1 1],'b')
hold on

% midpoint rule with Euler predictor
f = @(z) 1+z+z.^2/2;
Z = 1*(abs(f(X+1i*Y))<1);
subplot(3,2,2)
contourf(X,Y,Z,[1 1])
title('Midpoint w/ Euler Predictor','Interpreter','latex','FontSize',14)
xlabel('$Re(z)$','Interpreter','latex','FontSize',14)
ylabel('$Im(z)$','Interpreter','latex','FontSize',14)
grid on,xline(0),yline(0),axis square
subplot(3,2,6)
contour(X,Y,Z,[1 1],'r')

% Kutta's method
hk1 = @(z) z;
hk2 = @(z) z.*(1+hk1(z)/2);
hk3 = @(z) z.*(1-hk1(z)+2.*hk2(z));
f   = @(z) 1+hk1(z)/6 + 2*hk2(z)/3 + hk3(z)/6;
Z   = 1*(abs(f(X+1i*Y))<1);
subplot(3,2,3)
contourf(X,Y,Z,[1 1])
title("Kutta's method ",'Interpreter','latex','FontSize',14)
xlabel('$Re(z)$','Interpreter','latex','FontSize',14)
ylabel('$Im(z)$','Interpreter','latex','FontSize',14)
grid on,xline(0),yline(0),axis square
subplot(3,2,6)
contour(X,Y,Z,[1 1],'g')

% Heun's method
hk1 = @(z) z;
hk2 = @(z) z.*(1+hk1(z)/2);
hk3 = @(z) z.*(1+hk2(z)/2);
hk4 = @(z) z.*(1+hk3(z));
f   = @(z) 1 + hk1(z)/6 + hk2(z)/3 + hk3(z)/3 + hk4(z)/6;
Z   = 1*(abs(f(X+1i*Y))<1);
subplot(3,2,4)
contourf(X,Y,Z,[1 1])
title("Fourth order method ",'Interpreter','latex','FontSize',14)
xlabel('$Re(z)$','Interpreter','latex','FontSize',14)
ylabel('$Im(z)$','Interpreter','latex','FontSize',14)
grid on,xline(0),yline(0),axis square
subplot(3,2,6)
contour(X,Y,Z,[1 1],'c')

% DOPRI5
hk1 = @(z) z;
hk2 = @(z) z.*(1+1/5*hk1(z));
hk3 = @(z) z.*(1+3/40*hk1(z) + 9/40*hk2(z));
hk4 = @(z) z.*(1+44/45*hk1(z) - 56/15*hk2(z) + 32/9*hk3(z));
hk5 = @(z) z.*(1+19372/6561*hk1(z) - 25360/2187*hk2(z) + 64448/6561*hk3(z) - 212/729*hk4(z));
hk6 = @(z) z.*(1+9017/3168*hk1(z) - 355/33*hk2(z) + 46732/5247*hk3(z) + 49/176*hk4(z) - 5103/18656*hk5(z));
hk7 = @(z) z.*(1+35/384*hk1(z) + 500/1113*hk3(z) + 125/192*hk4(z) - 2187/6784*hk5(z) + 11/84*hk6(z));
subplot(3,2,5)
f = @(z) 1+5179/57600*hk1(z) + 7571/16695*hk3(z) + 393/640*hk4(z) - 92097/339200*hk5(z) + 187/2100*hk6(z) + 1/40*hk7(z);
Z   = 1*(abs(f(X+1i*Y))<1);
contourf(X,Y,Z,[1 1])
title("DOPRI5(4)",'Interpreter','latex','FontSize',14)
xlabel('$Re(z)$','Interpreter','latex','FontSize',14)
ylabel('$Im(z)$','Interpreter','latex','FontSize',14)
grid on,xline(0),yline(0),axis square
subplot(3,2,6)
contour(X,Y,Z,[1 1],'m')
xlabel('$Re(z)$','Interpreter','latex','FontSize',14)
ylabel('$Im(z)$','Interpreter','latex','FontSize',14)
grid on,xline(0),yline(0),axis square
legend('Forward Euler','Midpoint',"Kutta's",'Fourth order','DOPRI5(4)')