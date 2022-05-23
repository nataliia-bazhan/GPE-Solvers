function N = getN(obj, Psi)
N = sum(sum(sum(abs(Psi).^2))) * obj.grid.dV;
end