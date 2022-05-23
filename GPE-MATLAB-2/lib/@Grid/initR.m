function [r, rz, h, hz, dV] = initR(obj)
r = linspace(-obj.L/2, obj.L/2, obj.N);
rz = linspace(-obj.Lz/2, obj.Lz/2, obj.Nz);
h = r(2)-r(1);
hz = rz(2)-rz(1);
dV = h * h * hz;
end
