classdef Phase
   properties
      n = 3; % amount of sectors
      R_min = 0;  
      grid;
   end
   methods
      function obj = Phase(params)
          if nargin == 1
              if isfield(params, 'n')
                  obj.n = params.n;
              end
              if isfield(params, 'R_min')
                  obj.R_min = params.R_min;
              end
              if isfield(params, 'grid')
                  obj.grid = params.grid;
              end
          end
      end
   end
end