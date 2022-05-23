function [Psi, mu] = stepSGPE(obj, Psi, mu, t, projec)

xi_ra = randn(obj.grid.N, obj.grid.N, obj.grid.Nz)*obj.sigma; %noise generation
xi_rb = randn(obj.grid.N, obj.grid.N, obj.grid.Nz)*obj.sigma;
eta_r = (xi_ra + 1i*xi_rb)/...
        (obj.config.hbar*(1i - obj.config.gamma));
Psi = Psi + eta_r*obj.dt;
          
V = obj.getV(t);
          
Psi = ifftn(obj.ekk.*fftn(Psi));
Psi = exp((-1/obj.config.hbar) * obj.dt_f * ...
          (V + obj.config.g * abs(Psi.*conj(Psi)) - mu)).* Psi;
Psi = ifftn(obj.ekk.*fftn(Psi).*projec);

end