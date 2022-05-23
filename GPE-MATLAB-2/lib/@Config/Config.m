classdef Config
   properties
       hbar=1.054571800e-34; % Planck constant
       kb = 1.38064852e-23; % Boltzmann constant
       amu = 1.660539040e-27; % atomic mass unit
       a0=5.2917721067e-11; % Bohr radius
       
       td = 10; %sec, time of decay
       gamma = 0.05; %dissipation parameter
       T = 20e-9; %temperature
       N = 1.4*10^5;
       mu;
       M;
       as;
       g;
       
   end
   methods
      function obj = Config(params)
          if nargin == 1
              if isfield(params, 'N')
                  obj.N = params.N;
              end
              if isfield(params, 'ChE')
                  if params.ChE == 'Na23'
                      obj.M = 22.9897692809*obj.amu; % mass of sodium-23 atom
                      obj.as = 2.75e-9; % scattering length of sodium-23
                  elseif params.ChE == 'Rb87'
                      obj.M = 86.909180527*obj.amu; % mass of rhubidium-87 atom
                      obj.as = 100.87*obj.a0; %5.77e-9; % scattering length of rhubidium-87
                  end
              end
              if isfield(params, 'gamma')
                  obj.gamma = params.gamma;
              end
              if isfield(params, 'T')
                  obj.T = params.T;
              end
              if isfield(params, 'mu')
                  obj.mu = params.mu;
              end
              obj.g = (4*pi*obj.hbar^2*obj.as)/obj.M; %coupling strength
          end 
      end
   end
end