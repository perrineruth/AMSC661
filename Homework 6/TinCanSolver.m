% uxx + uyy = 0
close all
clc

% get a mesh of -pi <= x <= pi, 0 <= pi <= 2
nx = 60*10;
ny = 20*10;
hx = 2*pi/nx;
hy = 2/ny;
x = linspace(-pi,pi,nx+1);
y = linspace(2,0,ny+1); % easier for me to visualize y up top
[X,Y] = meshgrid(x,y);
% with left and bottom boundary removed
Xtrim = X(1:(end-1),2:(end));
Ytrim = Y(1:(end-1),2:(end));

% arrays
xvec = reshape(Xtrim,[],1);
yvec = reshape(Ytrim,[],1);
fvec = -max(cos(xvec),0);

% form matrix w/ blocks of size ny x ny
I = speye(ny);
T = spdiags([[ones(ny-2,1); 1; 0]/hy^2, -(2/hx^2+2/hy^2)*ones(ny,1), [0; 2; ones(ny-2,1)]/hy^2],-1:1,ny,ny);

IStructure = spdiags([ones(nx,1) zeros(nx,1) ones(nx,1)],-1:1,nx,nx);
IStructure(1,nx) = 1;
IStructure(nx,1) = 1;
A = kron(eye(nx),T) + kron(IStructure,I)/hx^2;

% solve
%T = timeit(@() A\fvec);
tic
u = A\fvec;
toc
umat = reshape(u,ny,nx);
ufull = [[umat(:,end), umat]; zeros(1,nx+1)];

% add in 

h1 = pcolor(X,Y,ufull);
set(h1, 'EdgeColor', 'none');
xlabel('x')
ylabel('y')
set(gca,'FontSize',18)
colorbar

figure()
h2 = surf(cos(X),sin(X),Y,ufull);
set(h2, 'EdgeColor', 'none');
axis equal
axis off
view([5 -5 5])

%clc
%T