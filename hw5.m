% ME 523 - HW #5
clc; clear; close all

solver_name = ["Gauss-Seidel","Successive Over-Relaxation","Conjugate Gradient"];
linestyle = ["-","--",":"];
% Create uniform mesh
nx=50;                         % Number of grid points in x
ny=50;                         % Number of grid points in y
x=linspace(0,1,nx+1);          % x-Grid (location of cell faces)
y=linspace(0,1,ny+1);          % y-Grid (location of cell faces)
xm=x(1:end-1)+(x(2)-x(1))/2;   % x-Grid (location of cell centers)
ym=y(1:end-1)+(y(2)-y(1))/2;   % y-Grid (location of cell centers)
dx=xm(2)-xm(1);                % Grid spacing in x
dy=ym(2)-ym(1);                % Grid spacing in y

% Generate RHS
b=zeros(nx*ny,1);
for i=1:nx
    for j=1:ny
        % Lexicographic ordering
        n=i+nx*(j-1);
        % Source (store at cell center)
        b(n)=exp(-((xm(i)-0.75)^2+(ym(j)-0.75)^2)/(0.05^2))...
            -exp(-((xm(i)-0.25)^2+(ym(j)-0.25)^2)/(0.05^2));
    end
end
b = b*dx*dy;

% Construct A
as = ones(nx*ny,1)*(dx/dy);
ae = ones(nx*ny,1)*(dy/dx);
aw = ones(nx*ny,1)*(dy/dx);
an = ones(nx*ny,1)*(dx/dy);
ap = -(as+ae+aw+an);

%check whether positive definite
B=[as an ap as an];
d=[nx*(ny-1) nx 0 -nx -(nx*(ny-1))];
A=(spdiags(B, d, nx*ny, nx*ny))';
for i = 1:nx*ny
    % Put ae in A
    if ~mod(i,nx)
        A(i-(nx-1),i)=ae(i);
    else
        A(i+1,i)=ae(i);
    end
    %Put aw in A
    if ~mod(i-1,nx)
        A(i+nx-1,i)=aw(i);
    else    
        A(i-1,i)=aw(i);
    end
end
K = full(A);
[~,p] = chol(A) % 0 if pd 

lim = 1e-7; % residual limit
maxIt = 10000; % iterational limit
omega = 1.5085; % relaxation factor
for j = 1:3
    res = [1];
    % Initialize phi
    phi_S=zeros(nx*ny,1);
    phi_0=zeros(nx*ny,1);
    phi_1=zeros(nx*ny,1);
    r = zeros(nx*ny,1);
    d = zeros(nx*ny,1);
    solver = j;
    if solver == 1
        % Gauss-Seidel
        while res(end)>lim  && length(res)<maxIt
            phi_0=phi_1;
            for i = 1:nx*ny
                phi_1(i)=-1/ap(i)*(as(i)*phi_1(1+mod(i-nx-1,ny*nx))...
                    + aw(i)*phi_1(i-1+(mod(i-1,nx)==0)*nx)...
                    + ae(i)*phi_1(i+1-(mod(i,nx)==0)*nx)...
                    + an(i)*phi_1(1+mod(i+nx-1,ny*nx)) + b(i));
            end
            res = [res norm(phi_1-phi_0)];
        end
    elseif solver == 2
        % SOR
        while res(end)>lim && length(res)<maxIt
            phi_0=phi_1;
            for i = 1:nx*ny
                phi_S(i)=-1/ap(i)*(as(i)*phi_S(1+mod(i-nx-1,ny*nx))...
                    + aw(i)*phi_S(i-1+(mod(i-1,nx)==0)*nx)...
                    + ae(i)*phi_1(i+1-(mod(i,nx)==0)*nx)...
                    + an(i)*phi_1(1+mod(i+nx-1,ny*nx)) + b(i));
            end
            phi_1 = phi_1 + omega*(phi_S-phi_1);
            res = [res norm(phi_1-phi_0)];
        end
    elseif solver == 3
        % Conjugate Gradient
        r = b - A*phi_1;
        rho = [norm(r)^2];
        while res(end)>lim && length(res)<maxIt
            phi_0=phi_1;
            if length(res)==1
                d = r;
            else
                beta = rho(end)/rho(end-1);
                d = r + beta*d;
            end
            e = A*d;
            a = rho(end)/(d'*e);
            phi_1 = phi_1 + a*d;
            r = r -a*e;
            rho = [rho r'*r];
            res = [res norm(phi_1-phi_0)];
        end
        phi_1=flip(phi_1);
    end
    
    figure(j)
    plot((nx*ny/2):1:((ny+2)*nx)/2,phi_1((nx*ny/2):1:((ny+2)*nx)/2));
    xlabel('x');
    ylabel('\phi');
    title(['Values of \phi at y=0.5 With ' solver_name(j)]);
    name = [solver_name(j)]
    print (name,'-dpdf','-fillpage')
    
    length(res)
    phi_1 = phi_1 - (1/(length(phi_1)))*sum(phi_1);
    sum(phi_1*dx*dy)
    figure(4)
    semilogy(1:length(res),res,linestyle(j),'LineWidth',2)
    title('Residuals');
    xlabel('Iterations');
    ylabel('\epsilon_k');
    hold on
end
xlim([0 1100])
hold off
legend('Gauss-Seidel','Successive Over-Relaxation','Conjugate Gradient')
print ('HW5_res','-dpdf','-fillpage')