function V = getV(obj, t)

V = zeros(obj.grid.N, obj.grid.N, obj.grid.Nz);
if ~isempty(obj.Vt)
    V = V + obj.Vt.getV();
end
if ~isempty(obj.Vb)
    V = V + obj.Vb.getV(t);
end
end