%% Initialization

units % load some common physical constants
omegar = 150*2*pi; % trap frequencies (150Hz and 600Hz)
omegaz = 600*2*pi;
tcoef = 1/(omegar); % time scale
rscale = sqrt(hbar/mRB/omegar); % spatial scale

N=128/4; % number of grid points for coordinate grid
xmax = 15e-6/rscale; %  (15 micrometers)  grid boundaries [-xmax,xmax]
x = gpuArray.linspace(-xmax,xmax,N);

Nz=32/4; % number of grid points for coordinate grid
zmax = (15e-6/rscale)/4;
z = gpuArray.linspace(-zmax,zmax,Nz);

grid = grid3d(xmax,N, xmax,N, zmax,Nz); 

V = @(X,Y,Z) 0.5*((X.^2 + Y.^2) + (omegaz/omegar)^2*Z.^2); % function representing the external potential (must have 3 arguments)
%task = ZNGtask(grid,V);
task = PGPEtask(grid, V);
%tasks = [ZNGtask(grid,V), PGPEtask(grid, V)]

Ntotal = 8.0e4;
task.Ntotal = Ntotal; % wave function normalization (total number of atoms)
task.g = 4*pi*aRB/rscale/sqrt(2*pi/omegaz*omegar); % nonlinear interaction constant
%task.T = 10;
%task.Vtd = @(o,t) 10*exp(-((o.grid.mesh.x+5*cos(t)).^2+(o.grid.mesh.y+5*sin(t)).^2));


%% Stationary state calculation

tstep = 0.02; % initial time step for imaginary time evolution
acc = 1e-6; % desired accuracy


for T = 1:1:165
    task.T = T;
    phi = task.groundstate_itp(tstep,acc); % phi will contain calculated stationary state
    
    cT = T*hbar*omegar/kb;
    cNc = task.grid.integrate(abs(task.init_state).^2);
    %cNt = task.grid_nt.integrate(task.init_state_nt);
    
    %imagesc(x, z, abs(squeeze(phi(:, N/2, :))).^2)
    %shg;
    %imagesc(x, z, transpose(squeeze(task.init_state_nt(:, N/2, :))))
    %colorbar;
    %shg;
    
    Tar(T) = cT;
    Nc(T) = cNc;
    %Nt(T) = cNt;
    
    display(['T = ', num2str(cT, '%.3g')])
    display([num2str(cNc)])
    %display([num2str(cNc), ', ', num2str(cNt), ', ', num2str(cNc + cNt)])
end

plot(Tar, Nc/Ntotal)
hold on

plot(Tar, Nt/Ntotal)
hold on

plot(Tar, (Nt + Nc)/Ntotal)
hold on 

legend('Nc', 'Nt', 'Ntot');
shg