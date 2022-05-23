function [phi, varargout] = geg2(task,dt,eps,phi0)

rscale = 8.8*10^(-7);

TT = task.T;
grid = task.grid;
grid_nt = task.grid_nt;
nx = grid.nx/2;
ny = grid.ny/2;
nz = grid.nz/2;
V = task.getVtotal(0);
g = task.g;
n_cn=task.n_crank;

NN0 = task.Ntotal;

if(nargin <= 3)
    phi0 = 'rand';
end

phi = sqrt(NN0)*grid.normalize(rand(size(grid.mesh.x),'like',V) + 1i*rand(size(grid.mesh.x),'like',V)); % random initial guess

nt = zeros(size(grid_nt.mesh.x),'like',V);
tmpext = task.vtrap_nt;%zeros(size(grid_nt.mesh.x),'like',V);
NNN = NN0;
NNt = 0;
ekk = exp(-grid.kk*dt);
MU = zeros(1000,1,'like',V);
MU2 = zeros(1000,1,'like',V);
delta = 1;
mu_old = 0;
i = 1;

tmp2 = real(phi.*conj(phi))*g+V;
while (delta > eps || mod(i,10)~=5)
    phi = exp(-tmp2*dt*0.5).*phi;
    phi = grid.ifft(ekk.*grid.fft(phi));
    phi = exp(-tmp2*dt*0.5).*phi;
    
	tmp = real(phi.*conj(phi));
    if(TT>0 && mod(i,10)==1 && task.Ntotal > 0) % for better performance and stability we do some initial iterations without a thermal cloud
        NNt = grid_nt.integrate(nt);
        NNN = NN0 - NNt;
        if(NNN<1)
            NNN=1;
            nt = nt*NN0/NNt; % we need to get the correct total number of particles even above Tc
        end
    end
    ncur = grid.integrate(tmp);
    mu = sqrt(NNN/ncur);
    MU(i) = mu;
    ncur = ncur*mu^2;

    phi=phi*mu;
	tmp = tmp*mu^2;
	tmp2 = tmp*g+V;
    MU2(i) = real(grid.integrate(conj(phi).*grid.lap(phi) + tmp.*(tmp2+2*g*nt(nx+1:3*nx,ny+1:3*ny,nz+1:3*nz))));
    HHk = real(grid.integrate(conj(phi).*grid.lap(phi)));
    HHi = real(grid.integrate(tmp.*(tmp2+2*g*nt(nx+1:3*nx,ny+1:3*ny,nz+1:3*nz))));
    %display(['Hk = ', num2str(HHk), ' Hi = ', num2str(HHi)])
    
    
    MU2(i) = MU2(i)/ncur;
    if(mod(i,10)==0)%i>150)
		if(TT>0)
		    tmpext(nx+1:3*nx,ny+1:3*ny,nz+1:3*nz) = V + 2*g*tmp;
            mmu = min(MU2(i),min(min(min(tmpext + 2*g*nt)))-1e-10); % compensate for possibly inaccurate chem. pot. calculation
            
		    ntt=(TT/(2*pi))^1.5*mypolylog(1.5,exp((mmu - tmpext - 2*g*nt)/TT)); % averaging increases stability for high temperatures
		    %display(['ssec = ', num2str(min(min(min(tmpext + 2*g*nt)))-1e-10)])
            %display(['mmu = ', num2str(mmu), ' mu2 = ', num2str(MU2(i))])%,' c3 = ', num2str(c3)])
            
            p = mypolylog(1.5,exp((mmu - tmpext - 2*g*nt)/TT));
            N=128/4; % number of grid points for coordinate grid
            xmax = 15e-6/(8.8*10^(-7)); %  (15 micrometers)  grid boundaries [-xmax,xmax]
            xxx = linspace(-xmax,xmax,N);
            imagesc(xxx, xxx, real(p(:, :, 4)))
    colorbar;
    drawnow;
            display(['Nt = ', num2str(NNt), '/', num2str(NNN)])
            if(NNN<=1)
		        NNt = grid_nt.integrate(ntt);
		        nt = (nt+ntt*NN0/NNt)*0.5;
		    else
		        nt=(nt+ntt)*0.5;
		    end
		end
        delta = abs(log(mu_old/mu))/dt^2/9;
        mu_old = MU(i-9);
    end
	tmp2 = tmp2 + 2*g*nt(nx+1:3*nx,ny+1:3*ny,nz+1:3*nz);
    i=i+1;
    if(i>=10000) 
        warning('Convergence not reached');
        break;
    end
end

task.init_state = phi;
task.init_state_nt = nt;

end
