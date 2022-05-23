classdef Grid
   properties
       GPU = 1; 
       L = 40*1.84*10^(-6);
       Lz = 10*1.84*10^(-6);
       N = 128*2;
       Nz = 64;
       
       r;
       rz;
       h;
       hz;
       dV;
       X;
       Y;
       Z;
       kk;
   end
   methods
      function obj = Grid(params)
          if nargin == 1
              if isfield(params, 'GPU')
                  obj.GPU = params.GPU;
              end
              if isfield(params, 'L')
                  obj.L = params.L;
              end
              if isfield(params, 'Lz')
                  obj.Lz = params.Lz;
              end
              if isfield(params, 'N')
                  obj.N = params.N;
              end
              if isfield(params, 'Nz')
                  obj.Nz = params.Nz;
              end
          end
          [obj.r, obj.rz, obj.h, obj.hz, obj.dV] = obj.initR();
          [obj.X, obj.Y, obj.Z] = obj.initXYZ();
          obj.kk = obj.initKK();
      end
   end
end