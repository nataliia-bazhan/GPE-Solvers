function [X, Y, Z] = initXYZ(obj)
[X, Y, Z] = meshgrid(obj.r, obj.r, obj.rz);
if obj.GPU
    X = gpuArray(X);
    Y = gpuArray(Y);
    Z = gpuArray(Z);
end
end

