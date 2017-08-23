%parameter.plot_onoff = 'off';  % turn on plot when you don't run PaperFigure.m

%% Parameter Initialization
% You have to tune the tradeoff parameter lambda!!
parameter.prop  = 1; 
parameter.MAX_ITER_Reweight =5;
parameter.Gammai_threshold = 1*[1e-3  1e-3  1e-3 1e-3 1e-2 1e-2 1e-2 1e-1 1e-1];
parameter.diff_onoff = 'off';
parameter.truncate = 'on';
parameter.tmp = 'off';
parameter.info = cpuinfo;
parameter.core = 5;
parameter.rho_update = 'off';
parameter.rho =  1;
parameter.alpha = 0.7;
parameter.epsilon = 1e-10;
parameter.scale = 1e0;
parameter.scale_w = 1e0;
parameter.dimblock = 5e3; %dimension of the each block
parameter.central_threshold = 5e3;   % threshold for using centralised algorithm even distributed_onoff is 'on'
parameter.Anorm_onoff = 'on';
diaginv = {'woodbury', 'lanzcos', 'full'};
prune = {'gamma','sensitivity','rank'};
method = {'woodbury', 'cholesky', 'lcg', 'full'};
parameter.diaginv =diaginv{1};
parameter.prune = prune{1};
parameter.method = method{1};
parameter.group_type = 'SMV';
parameter.multiple = 'on';
parameter.single = 'off';
parameter.group_onoff = parameter.multiple;
parameter.derivative ='solver'; % change this accordingly
parameter.derivativekth = 1;

prop = parameter.prop;
distributed_onoff = 'off'; % If you switch to admm solver, you can enable distributed computation by specifying multi cores
solver = {'cvx','admm'}; % Switch which solver you want to use; 
                                             % CVX is using CVX toolbox calling Sedumi;
                                             % admm is customed implemented, which is faster
solver = solver{1};
parameter.solver = solver;