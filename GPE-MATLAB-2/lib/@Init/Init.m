classdef Init
   properties
       config; 
       grid;
       Vt; 
       Vb;
       
       phase;
       
       dt_itp = 2*10^(-5); % time step for ITP
       MU_DIF = 10^(-33); % criterion of breaking the cycle
       n_iter = 500;
       t = 0;
       
       ekk;
   end
   methods
      function obj = Init(params)
          if nargin == 1
              if isfield(params, 'grid')
                  obj.grid = params.grid;
              end
              if isfield(params, 'config')
                  obj.config = params.config;
              end
              if isfield(params, 'phase')
                  obj.phase = params.phase;
              end
              if isfield(params, 'Vt')
                  obj.Vt = params.Vt;
              end
              if isfield(params, 'Vb')
                  obj.Vb = params.Vb;
              end
              if isfield(params, 'n_iter')
                  obj.n_iter = params.n_iter;
              end
          end 
          obj.ekk = exp(-0.5 * obj.grid.kk * obj.dt_itp * ...
                       (obj.config.hbar / obj.config.M));
      end
      
      
   end
end