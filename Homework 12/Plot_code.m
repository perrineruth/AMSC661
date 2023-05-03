%% plot 1.b)
% Perrin Ruth
close all
clc
fig = figure(); 

f = @(p) -p.*log(p);
df = @(p) -log(p)-1;
s = (f(0.1)-f(0.9))/(0.1-0.9);


t = [0,1];
plot([0 s],t,'k--','LineWidth',2)
hold on

for i = 1:40
   % left characteristic
   if -i/5 + df(.1) > s
       tintersect = -i/5/(s-df(.1));
       plot([-i/5, -i/5+df(.1)*tintersect],[0,tintersect],'k')
   else
       plot([-i/5, -i/5+df(.1)],t,'k')
   end
    
   if i/10 + df(.9) < s
       tintersect = i/10/(s-df(.9));
       plot([i/10, i/10+df(.9)*tintersect],[0,tintersect],'k')
   else
       plot([i/10, i/10+df(.9)],t,'k')
   end
end


xlabel('x')
ylabel('t')
fontsize(fig, 12,'points')
xlim([-1/2 1/2])
ylim([0 1])

%% Part 1.c
close all
hold on
xs = fsolve(@(x)1+(5*pi/9+atan(x)).*2*x,10);
Tb = (5*pi/9+atan(xs))*(1+xs^2)
p0 = @(x) 0.5+0.9/pi*atan(x);
for x0 = linspace(-10,10,200)
    slope = df(p0(x0));
    plot([x0,x0+slope*2],[0, 2],'k')
end
plot([-1 1],[Tb Tb],'k--','LineWidth',2)
xlim([-1 1])

s = (f(.05)-f(.95))/(.05-.95)