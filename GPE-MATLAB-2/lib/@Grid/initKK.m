function kk = initKK(obj)
k = [(0 : obj.N/2 - 1)*2*pi/obj.L -(obj.N/2 : -1 : 1)*2*pi/obj.L]; 
kz = [(0 : obj.Nz/2 - 1)*2*pi/obj.Lz -(obj.Nz/2 : -1 : 1)*2*pi/obj.Lz];
[KX,KY,KZ] = meshgrid(k,k,kz);
          
if obj.GPU
    kk = gpuArray((KX.^2+KY.^2+KZ.^2)/2); 
else
    kk = ((KX.^2+KY.^2+KZ.^2)/2);
end
end