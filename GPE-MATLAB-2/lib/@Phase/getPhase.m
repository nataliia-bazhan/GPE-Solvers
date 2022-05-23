function phaseFilter = getPhase(obj)
phaseFilter = ones(obj.grid.N, obj.grid.N, obj.grid.Nz);

for ii = 1:obj.n
    rand_phase = 2*pi*rand;
    filter = ((atan2(obj.grid.Y, obj.grid.X) > ...
                                 -pi + ii*2*pi/obj.n) & ...
              (atan2(obj.grid.Y, obj.grid.X) < ...
                                 -pi + (ii + 1)*2*pi/obj.n)).*...
              ((obj.grid.X.^2 + obj.grid.Y.^2) > obj.R_min.^2);
          
    phaseFilter = phaseFilter.*exp(1i*rand_phase*filter);
end
end