function KdVrkm()
close all
fsz = 20; % fontsize
% solves u_t + u_{xxxx} + u_{xx} + (0.5u^2)_x = 0, i.e.,
% u_t = -u_{xxxx} - u_{xx} - (0.5u^2)_x


N = 512;
L = 32*pi;
x = linspace(0,32*pi,N+1);
x(N + 1) = [];
k = -N/2 : (N/2 - 1); % wave numbers
u = zeros(1,N);
% initial data
u0 = cos(x/16).*(1+sin(x/16));

dt = 0.1; % time step
figure; clf; 
hpic = plot(x,u0,'LineWidth',2,'color','r'); % plot of the numerical solution
hold on;
grid
xlim([0,L]);
set(gca,'Fontsize',fsz);
xlabel('x','FontSize',fsz);
ylabel('u','FontSize',fsz);
drawnow
%
tmax = 200;
t = 0;
freq = k.*(2*pi/L); % frequencies
freq2 = freq.^2;
freq4 = freq.^4;
forwardexp=exp((freq2-freq4)*dt); % in the Fourier space, uhat = e3.*vhat
tvec = 0:dt:200;
UData = zeros(length(x),length(tvec));
UData(:,1) = u0;
index = 1;
while (t<tmax) 
    index = index + 1;
    t=t+dt;
    vhat=fftshift(fft(u0)); % v in the Fourier space
    % RK4 step in the Fourier space
    k1=rhs(0,vhat);
    k2=rhs(0.5*dt,vhat+0.5*dt*k1);
    k3=rhs(0.5*dt,vhat+0.5*dt*k2);
    k4=rhs(dt,vhat+dt*k3);
    vhat_new=vhat+dt*(k1+2*k2+2*k3+k4)/6;
    % return to the original space and the original variable u
    unew=ifft(ifftshift(forwardexp.*vhat_new)); % return to u in the x-space
    set(hpic,'xdata',x,'ydata',real(unew));
    u0=unew;
    drawnow
    UData(:,index) = real(u0);
end
close all
imagesc(tvec,x,UData)
xlabel('time (t)')
ylabel('space (x)')
end
%%
function RHSvhat=rhs(dt,vhat)
% v should be a row vector
% RHSvhat = - e^{-tL}(1i*k*hat{(e^{tL}v)^2/2} 
N=size(vhat,2);
L = 32*pi;
k=-N/2 : (N/2 - 1);
freq =k.*(2*pi/N);
freq2 = freq.^2;
freq4 = freq.^4;
forwardexp=exp((freq2-freq4)*dt);
backwardexp=exp(-(freq2-freq4)*dt);
vhat1=vhat.*forwardexp;          % e^{tL}v in the Fourier space 
v1=ifft(ifftshift(vhat1));      % exp(tL)v in the x-space
v2=0.5*v1.^2;          % [exp(tL)v]^2 in the x-space
RHSvhat=-backwardexp.*(1i*freq).*fftshift(fft(v2)); % exp(-tL)[[(exp(tL)v)]_x] in the Fourier space
end
