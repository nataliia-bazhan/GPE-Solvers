function [Psi, mu] = getITP(obj)
rng('shuffle');
%tic
phi0 = rand(obj.grid.N, obj.grid.N, obj.grid.Nz) + ...
            1i*rand(obj.grid.N, obj.grid.N, obj.grid.Nz);
             
cons = sqrt(obj.config.N/(sum(sum(sum(abs(phi0).^2)))*obj.grid.dV));
Psi = phi0.*cons;
    
MU = zeros(obj.n_iter, 1);
NN_ar = zeros(obj.n_iter, 1);

V = obj.getV(obj.t);

for i=1:obj.n_iter
    NN_ar(i) = i;
              
    Psi = ifftn(obj.ekk.*fftn(Psi));
    Psi = exp(-(obj.dt_itp/(obj.config.hbar))*...
                  (V + obj.config.g*abs(Psi.*conj(Psi)))).*Psi;
    Psi = ifftn(obj.ekk.*fftn(Psi));
    Psi2 = abs(Psi.*conj(Psi));
    
    imagesc(obj.grid.r, obj.grid.r, abs(Psi(:, :, 1)))
    hold on;
    drawnow;
              
    cons = sqrt(obj.config.N/(sum(sum(sum(Psi2)))* obj.grid.dV)); 
    Psi = Psi*cons;
    Psi2 = Psi2*cons^2;
    MU(i) = obj.config.hbar/obj.dt_itp*log(gather(cons));
              
    if i > 1
        dMU = abs(gather(MU(i) - MU(i-1)))/obj.dt_itp;
        if dMU < obj.MU_DIF
            last_ind = i;
            break;
        end
    end
    last_ind = i;
end

if ~isempty(obj.phase)
    Psi = Psi.*obj.phase.getPhase();
end

%% Residual
Hphi = ifftn((obj.config.hbar^2/obj.config.M)*obj.grid.kk.*fftn(Psi)) + ...
             (V + obj.config.g*Psi2).*Psi;
mu = gather(sum(sum(sum(abs(conj(Psi).*Hphi).*obj.grid.dV))))./obj.config.N;

display(['ITP finished with mu = ', num2str(mu)])
end