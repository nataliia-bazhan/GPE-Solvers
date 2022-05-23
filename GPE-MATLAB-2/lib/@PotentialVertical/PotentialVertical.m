classdef PotentialVertical
   properties
      w = 0;
      n = 3;
      R_min = 0;
      U_max = 10.6*2*pi*123*1.054*10^(-34);
      a_max = 0.6*1.89*10^-6;
      t1 = 0.01;
      t2 = 0.01;
      grid;
      config; 
   end
   methods
      function obj = PotentialVertical(params)
          if nargin == 1
              if isfield(params, 'grid')
                  obj.grid = params.grid;
              end
              if isfield(params, 'config')
                  obj.config = params.config;
              end
              if isfield(params, 't1')
                  obj.t1 = params.t1;
              end
              if isfield(params, 't2')
                  obj.t2 = params.t2;
              end
              if isfield(params, 'w')
                  obj.w = params.w;
              end
              if isfield(params, 'n')
                  obj.n = params.n;
              end
              if isfield(params, 'R_min')
                  obj.R_min = params.R_min;
              end
              if isfield(params, 'a_max')
                  obj.a_max = params.a_max;
              end
              if isfield(params, 'U_max')
                  obj.U_max = params.U_max;
              end
          end
      end
      
      
      
      
      
   end
end