%% Problem 1b - ODE solvers on the Van Der Pol Eq.
clf
% Equation
f = @(t,y,u) [y(2); u*((1-y(1).^2).*y(2))-y(1)];

tmax = 1000;  % max time for sim
y0   = [2;0]; % init conditions

% get 3x2 plots (col gives solver)
for i = 1:3
    u = 10^i; % VDP const
    % plot solutions
    % Runge Kutta
    [t,y] = ode15s(@(t,y) f(t,y,u),[0,tmax],y0);
    subplot(3,3,3*i-2);
    plot(y(:,1),y(:,2))
    ylabel('$y_2$','Interpreter','latex','FontSize',20)
    if i == 1
        title('Solution','Interpreter','latex','FontSize',20)
    end
    if i == 3
        xlabel('$y_1$','Interpreter','latex','FontSize',20)
    end

    % compute times
    eps_vec = [1e-6,1e-9,1e-12];
    time_vec_RK    = zeros(1,3);
    time_vec_stiff = zeros(1,3);
    for k = 1:3
        eps = eps_vec(k);
        % with Runge Kutta
        options = odeset('RelTol',eps,'AbsTol',eps);
        tstart = cputime;
        [~,~] = ode45(@(t,y) f(t,y,u),[0 tmax],y0,options);
        tend = cputime;
        time_vec_RK(k) = tend - tstart;
        tstart = cputime;
        [~,~] = ode15s(@(t,y) f(t,y,u),[0 tmax],y0,options);
        tend = cputime;
        time_vec_stiff(k) = tend - tstart;
    end
    % plot RK times
    subplot(3,3,3*i-1);
    loglog(eps_vec,time_vec_RK)
    ylabel('time (T)','Interpreter','latex','FontSize',20)
    if i == 1
        title('Runge Kutta (ode45)','Interpreter','latex','FontSize',20)
    end
    if i == 3
        xlabel('Error ($\epsilon$)','Interpreter','latex','FontSize',20)
    end
    % plot stiff times
    subplot(3,3,3*i);
    loglog(eps_vec,time_vec_stiff)
    ylabel('time (T)','Interpreter','latex','FontSize',20)
    if i ==1
        title('Stiff (ode15s)','Interpreter','latex','FontSize',20)
    end
    if i == 3
        xlabel('Error ($\epsilon$)','Interpreter','latex','FontSize',20)
    end
end

%% Problem 1c - Arenstorf equation
clf
% params
u = 0.012277471;
y0 = [0.994; 0; 0; -2.001585106];
f = @(t,y) diff_Arenstorf(t,y,u);
Tmax =  17.0652165601579625588917206249;
eps = 1e-12;
options = odeset('RelTol',eps,'AbsTol',eps);

% plot single orbit
subplot(3,1,1)
[~,y] = ode45(f,[0,Tmax],y0,options);
plot(y(:,1),y(:,2))
title('Single cycle $T_{\max}$ = 17.07','interpreter','latex','FontSize',14)
axis square

% for Tmax = 100
subplot(3,1,2)
[~,y] = ode45(f,[0,100],y0,options);
plot(y(:,1),y(:,2))
title('ode45 $T_{\max}$=100','interpreter','latex','FontSize',14)
axis square

subplot(3,1,3)
[~,y1] = ode23s(f,[0,100],y0,options);
[~,y2] = ode89(f,[0,100],y0);
[~,y3] = ode45(f,[0,100],y0);
plot(y1(:,1),y1(:,2),y2(:,1),y2(:,2),y3(:,1),y3(:,2))
title('varying methods $T_{\max}$=100','interpreter','latex','FontSize',14)
lgd = legend('ode23s','ode89','ode45');
lgd.FontSize = 14;
xlim([-3 3])
ylim([-3 3])
axis square


% diff eq
function [Dy] = diff_Arenstorf(t,y,u)
    r1 = ((y(1)+u)^2+y(2)^2)^(1/2);
    r2 = ((y(1)-1+u)^2+y(2)^2)^(1/2);
    Dy = [y(3);
        y(4);
        y(1) + 2*y(4) - (1-u)*(y(1)+u)/r1^3 - u*(y(1)-1+u)/r2^3;
        y(2) - 2*y(3) - (1-u)*y(2)/r1^3 - u*y(2)/r2^3];
end
