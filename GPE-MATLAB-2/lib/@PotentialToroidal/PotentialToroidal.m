classdef PotentialToroidal
   properties
      wr = 2*pi*0.5*10^(3);
      wz = 2*pi*1.58*10^(3); % z trap frequence
      R = 19.23*10^(-6);  
      R2;
      nRings = 2; %amount of the coaxial rings
      grid;
      config; 
   end
   methods
      function obj = PotentialToroidal(params)
          if nargin == 1
              if isfield(params, 'wr')
                  obj.wr = params.wr;
              end
              if isfield(params, 'wz')
                  obj.wz = params.wz;
              end
              if isfield(params, 'R')
                  obj.R = params.R;
              end
              if isfield(params, 'R2')
                  obj.R2 = params.R2;
              end
              if isfield(params, 'nRings')
                  obj.nRings = params.nRings;
              end
              if isfield(params, 'grid')
                  obj.grid = params.grid;
              end
              if isfield(params, 'config')
                  obj.config = params.config;
              end
          end
      end
   end
end