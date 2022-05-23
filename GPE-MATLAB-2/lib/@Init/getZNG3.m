function [Psi, nt] = getZNG3(obj)

V = obj.getV(obj.t);

% initial Phi
if isempty(obj.config.mu)
    phi0 = rand(obj.grid.N, obj.grid.N, obj.grid.Nz) + ...
                1i*rand(obj.grid.N, obj.grid.N, obj.grid.Nz);  
    Psi = sqrt(obj.config.N/(sum(sum(sum(abs(phi0).^2)))*...
           obj.grid.dV))*phi0;
else
    Psi = real(sqrt(complex(obj.config.mu - V)/obj.config.g)); % use only Thomas-Fermi approximation as initial guess if mu_init is set
end

nt = zeros('like', Psi);

% initial amount of thermal and condensed atoms
if isempty(obj.config.mu)
    Nc = obj.config.N;
else
    Nc = 1;
end
Nt = 0;

NT = zeros(100, 1);
C = zeros(1000, 1);
delta = 1;
eps = 10^(-4);
i = 1;

Nl = V + obj.config.g*real(Psi.*conj(Psi)) + 2*obj.config.g*nt;

while (delta > eps)
    Psi = exp(-(obj.dt_itp/(2*obj.config.hbar))*Nl).*Psi;
    Psi = ifftn(obj.ekk.*fftn(Psi));
    Nl = V + obj.config.g*abs(Psi.*conj(Psi)) + 2*obj.config.g*nt;
    Psi = exp(-(obj.dt_itp/(2*obj.config.hbar))*Nl).*Psi;
   
    if(obj.config.T > 0 && mod(i,10) == 1 && isempty(obj.config.mu)) % for better performance and stability we do some initial iterations without a thermal cloud
        Nt = sum(sum(sum(nt.*obj.grid.dV)));
        Nc = obj.config.N - Nt;
        if(Nc < 1)
            Nc = 1;
            nt = nt*obj.config.N/Nt; % we need to get the correct total number of particles even above Tc
        end
    end
    
    N = sum(sum(sum(real(Psi.*conj(Psi)).*obj.grid.dV)));
    
    if isempty(obj.config.mu)
        c = sqrt(Nc/N);
    else
        c = exp((obj.config.mu/obj.config.hbar)*obj.dt_itp);
    end
    N = N*c^2;
    C(i) = gather(c);

    Psi=Psi*c;
	Nl = V + obj.config.g*real(Psi.*conj(Psi)) + 2*obj.config.g*nt;
    
    Hk = conj(Psi).*(obj.config.hbar^2/(2*obj.config.M)).*...
               ifftn(obj.grid.kk.*fftn(Psi));
    Hi = Psi.*conj(Psi).*Nl;
    mu2 = real(sum(sum(sum((Hk + Hi).*obj.grid.dV))))/N;

    MU2(i) = gather(mu2);
    
    if(mod(i,10) == 0)
		if(obj.config.T > 0)
            mmu = min(MU2(i),min(min(min(Nl)))); % compensate for possibly inaccurate chem. pot. calculation
            
            ntt = (obj.config.kb*obj.config.T*obj.config.M/(2*pi))^(3/2)*...
                  (obj.grid.dV/obj.config.hbar^3)*...
                   myPolylog(1.5, exp((mmu - Nl)/(obj.config.kb*obj.config.T)));
               
            mmu
            %imagesc(obj.grid.r, obj.grid.r, real(p(:, :, 4)))
            %colorbar;
            %drawnow;
            
            Nt = sum(sum(sum(ntt)));
            NT(i/10) = gather(Nt);
            display(['Nt = ', num2str(Nt), '/', num2str(obj.config.N)])
            
            if(Nc <= 1 && isempty(obj.config.mu))
		        nt = (nt+ntt*obj.config.N/Nt)*0.5;
		    else
		        nt=(nt+ntt)*0.5;
		    end
        end
        if i > 19
            delta = abs(NT(i/10) - NT((i/10) - 1));
        end

    end
	Nl = V + obj.config.g*real(Psi.*conj(Psi)) + 2*obj.config.g*nt;
    i=i+1;
    if(i>=10000) 
        warning('Convergence not reached');
        break;
    end
end


end
