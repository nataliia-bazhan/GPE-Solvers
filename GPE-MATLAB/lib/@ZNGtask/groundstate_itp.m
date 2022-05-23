function [phi, varargout] = geg(task,dt,eps,phi0)

TT = task.T;
grid = task.grid;
grid_nt = task.grid_nt;
nx = grid.nx/2;
ny = grid.ny/2;
nz = grid.nz/2;
V = task.getVtotal(0);
g = task.g;
omega = task.omega;
n_cn=task.n_crank;

NN0 = task.Ntotal;

if(nargin <= 3)
    phi0 = 'rand';
end

phi = sqrt(NN0)*grid.normalize(rand(size(grid.mesh.x),'like',V) + 1i*rand(size(grid.mesh.x),'like',V)); % random initial guess

nt = zeros('like',V);
tmpext = task.vtrap_nt;%zeros(size(grid_nt.mesh.x),'like',V);
NNN = NN0;
NNt = 0;
ekk = exp(-grid.kk*dt);
MU = zeros(1000,1,'like',V);
MU2 = zeros(1000,1,'like',V);
delta = 1;
mu_old = 0;
i = 1;

    size(real(phi.*conj(phi)))
    size(V)
    size(nt)

tmp2 = real(phi.*conj(phi))*g+V + 2*g*nt;

while (delta > eps || mod(i,10)~=5)
    phi = exp(-tmp2*dt*0.5).*phi;
    phi = grid.ifft(ekk.*grid.fft(phi));
    phi = exp(-tmp2*dt*0.5).*phi;
    
	tmp = real(phi.*conj(phi));
    if(TT>0 && mod(i,10)==1) % for better performance and stability we do some initial iterations without a thermal cloud
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
	tmp2 = real(phi.*conj(phi))*g + V + 2*g*nt;
    MU2(i) = real(grid.integrate(conj(phi).*grid.lap(phi) + ...
             tmp.*(tmp2)))/ncur;
    i
    MU2(i)
    if(mod(i,10)==0)%i>150)
		if(TT>0)
		    %tmpext(nx+1:3*nx,ny+1:3*ny,nz+1:3*nz) = V + 2*g*tmp;
            mmu = min(MU2(i),min(min(min(tmp2)))-1e-10); % compensate for possibly inaccurate chem. pot. calculation
            display('.........')
            mmu
		    ntt=(TT/(2*pi))^1.5*mypolylog(1.5,exp((mmu - tmp2)/TT)); % averaging increases stability for high temperatures
		    if(NNN<=1)
		        NNt = grid_nt.integrate(ntt);
		        nt = (nt+ntt*NN0/NNt)*0.5;
            else
                NNt
		        nt=(nt+ntt)*0.5;
                display('+++++')
		    end
		end
        delta = abs(log(mu_old/mu))/dt^2/9;
        mu_old = MU(i-9);
    end
	tmp2 = real(phi.*conj(phi))*g+V + 2*g*nt;
    i=i+1;
    if(i>=10000) 
        warning('Convergence not reached');
        break;
    end
end
i=i-1;
if(nargout >= 2)
    MU = real(MU(1:nnz(MU)));
    if(task.Ntotal > 0)
        MU = 1/dt * log(MU);
        task.current_mu = MU(end);
        task.current_n = task.Ntotal;
    end
    varargout{1} = MU;
end
if(nargout >= 3)
    MU2 = real(MU2(1:nnz(MU2)));
%     MUEX = MU2(i) - (MU2(i)-MU2(i-1))^2/(MU2(i)-2*MU2(i-1)+MU2(i-2)); % exponential extrapolation
    varargout{2} = MU2;
end

task.init_state = phi;
task.init_state_nt = nt;

end
