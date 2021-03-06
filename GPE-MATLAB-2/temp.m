clear;
config = Config(struct('N', 6*10^5,...
                       'ChE', 'Na23'));
               
grid = Grid(struct('N', 32,...
                   'Nz', 8));
              
Vt = PotentialToroidal(struct('wr', 2*pi*123, ...
                              'wz', 2*pi*600, ...
                              'grid', grid, ...
                              'config', config));
                          
Vb = PotentialVertical(struct('grid', grid, ...
                              'config', config));


[Psi, mu] = Init(struct('grid', grid, ...
                        'config', config, ...
                        'Vt', Vt, ...
                        'Vb', Vb)).getITP(); 
                 
                   
dyn = Dynamics(struct('grid', grid, ...
                      'config', config, ...
                      'Folder', 'new folder',...
                      'Vt', Vt, ...
                      'Vb', Vb));
                  
t = 0;
dyn.runGPE(Psi, mu, t)








