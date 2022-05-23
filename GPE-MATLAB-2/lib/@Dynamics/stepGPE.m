function [Psi, mu] = stepGPE(obj, Psi, mu, t)

V = obj.getV(t);
N_prev = obj.getN(Psi);

Psi = ifftn(obj.ekk.*fftn(Psi));
Psi = exp((-1/obj.config.hbar) * obj.dt_f * ...
          (V + obj.config.g * abs(Psi.*conj(Psi)) - mu)).* Psi;
Psi = ifftn(obj.ekk.*fftn(Psi));
          
Nbase = obj.getN(Psi);
G = -(1/(obj.dt)) * log(gather (Nbase)/gather(N_prev));
mu = mu + obj.config.hbar*(G - obj.inv_td) * ...
    (1 + obj.config.gamma^2)/(2 * obj.config.gamma);

end