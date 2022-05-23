function saveData(obj, Psi, s)
ts = obj.ts;
save([obj.Data 't.mat'], 'ts');
          
ls = obj.ls;
save([obj.Data 'l.mat'], 'ls');
          
mus = obj.mus;
save([obj.Data, 'mu.mat'], 'mus');

save([obj.PsiMat 'Psi' num2str(s,'%.0f') '.mat'],'Psi');
end