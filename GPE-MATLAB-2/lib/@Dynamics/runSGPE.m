function runSGPE(obj, Psi, mu, t)
ecut = mu + obj.config.kb*obj.config.T*log(2);
kcut = sqrt(2*obj.config.M*ecut/obj.config.hbar^2);
projec = obj.grid.kk*2<=kcut^2;
          
s = 0;
for kj = 0 : obj.nit
    if (kj == obj.n_count_integrals*s)
        s = s + 1;

        obj.ts(s) = t;
        obj.ls(s) = obj.getL(Psi);
        obj.mus(s) = mu;
                  
        obj.saveData(Psi, s);
        obj.savePic(Psi, s)
        display(['t = ' num2str(t,'%.3f')]);
    end
    [Psi, mu] = obj.stepSGPE(Psi, mu, t, projec);
    t = t + obj.dt;
end
end