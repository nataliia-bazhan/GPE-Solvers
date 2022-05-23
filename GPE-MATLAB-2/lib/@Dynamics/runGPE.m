function runGPE(obj, Psi, mu, t)

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
    
    [Psi, mu] = obj.stepGPE(Psi, mu, t);
    t = t + obj.dt;
end

end
      