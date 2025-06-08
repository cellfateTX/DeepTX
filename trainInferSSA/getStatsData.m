root_path = "/public/home/zhjiajun/hzw/academic_code/DeepGTM/deepGTM/data/simulationData";
parmas_path = fullfile(root_path,"params_GTM.mat");
saved_path = "stats_params_GTM.mat";

parmas_ma = load(parmas_path);
params = parmas_ma.sim_paramsGTM;
[m,n]=size(params);
stat_the = [1,1,1,1,1,1];
for i=1:n
param_true.kon = params(3,i);
param_true.ron = params(4,i);
param_true.koff = params(1,i);
param_true.roff = params(2,i);
param_true.mu = params(5,i);
param_true.delta = 1;
param_true.x0 = [1,0,0];
param_true.tottime = 200;


% Simulation algorithm for GTM.
statis_ther = statisGTM(param_true,4);
stat_the = [stat_the;statis_ther];
sprintf('stat_%d.mat',i); 
end
stat_the = stat_the(2:end,:);
file_name = fullfile(root_path, saved_path);
file_name
save(file_name,'stat_the','params');
"done"

