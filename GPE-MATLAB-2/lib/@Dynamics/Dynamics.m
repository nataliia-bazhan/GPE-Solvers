classdef Dynamics
   properties
      grid;
      config; 
      
      Folder;
      Pictures;
      Data;
      PsiMat;
      
      dt = 5*10^-6;
      n_count_integrals = 200;
      t_init = 0.0;
      t_fin = 0.251;
      
      nit;
      dt_f;% time step (with gamma included)
      inv_td;
      ekk;
      
      sigma;
      
      Vt;
      Vb;
      
      ts;
      ls;
      mus;

   end
   methods
      function obj = Dynamics(params)
          if nargin == 1
              if isfield(params, 'grid')
                  obj.grid = params.grid;
              end
              if isfield(params, 'config')
                  obj.config = params.config;
              end
              if isfield(params, 'Vt')
                  obj.Vt = params.Vt;
              end
              if isfield(params, 'Vb')
                  obj.Vb = params.Vb;
              end
              if isfield(params, 'Folder')
                  obj.Folder = [params.Folder, '/'];
              end
          end
          
          obj.dt_f = obj.dt*1i/(1+1i*obj.config.gamma);% time step (with gamma included)
          obj.ekk = exp(-0.5 *(obj.config.hbar/obj.config.M) * ...
                        obj.grid.kk*obj.dt_f);
          obj.inv_td = 1/(obj.config.td);
          obj.nit= floor((obj.t_fin - obj.t_init)/obj.dt);
          
          obj.ts = zeros(1, floor(obj.nit/obj.n_count_integrals));
          obj.ls = zeros(1, floor(obj.nit/obj.n_count_integrals));
          obj.mus = zeros(1, floor(obj.nit/obj.n_count_integrals));
          
          obj.sigma = sqrt(obj.config.gamma*...
                          (obj.config.hbar*obj.config.kb*obj.config.T)/...
                          (obj.grid.dV*obj.dt));
         
          obj.Pictures = [obj.Folder, 'Pictures/'];
          obj.Data = [obj.Folder, 'Data/'];
          obj.PsiMat = [obj.Folder, 'Psi/'];
          
          if ~exist(obj.Pictures, 'dir')
              mkdir(obj.Pictures)
          end
          if ~exist(obj.Data, 'dir')
              mkdir(obj.Data)
          end
          if ~exist(obj.PsiMat, 'dir')
              mkdir(obj.PsiMat)
          end
      end
   end
end