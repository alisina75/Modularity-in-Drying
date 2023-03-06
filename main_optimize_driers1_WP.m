clc;clear;
Xwb = 0.75; %Input Moisture
Tmax = 70; Tmin = 30;
time_l = 0.01; time_u = 50;
K_stages= 5;
Mwb_des = 0.07;
num_paths = 2^K_stages; flag = 0;
paths = dec2bin(0:2^K_stages-1) - '0'; %size : 2^K_stages * K_stages 
pert = 0.0000000001;

% Setting the DA parameters
bmax=100; alpha=1.15; b = 0.01;
counter = ceil(log(bmax/b)/log(alpha));
res_time = time_l * ones(size(paths)); 
T = Tmin * ones(size(paths));
V = (1/num_paths)*ones(1,num_paths); max_prob = max(V);
t0 = reshape(res_time',[1 numel(paths)]); % size: 1 * (K_stages * 2^K_stages)
T0 = reshape(T',[1 numel(paths)]); % size: 1 * (K_stages * 2^K_stages)
x0 = [t0, T0]'; % size: (2 * K_stages * 2^K_stages) *1

%% annealing part
while max_prob < 0.999
    disp(['b = ',num2str(b)]);
    disp(['counter = ',num2str(counter)])
    s = numel(x0)/2;
    fun = @(x)free_energy1(x,b,V,Xwb);
    nonlcon = @(x)nonlin_const1(x,V,Xwb);
    lbtime = time_l * ones(s,1);
    lbT = Tmin * ones(s,1);
    lb = [lbtime;lbT];
    ubtime = time_u * ones(s,1);
    ubT = Tmax * ones(s,1);
    ub = [ubtime;ubT];
    options = optimoptions('fmincon','MaxFunctionEvaluations',3*1e04,...
    'ScaleProblem' , true, ...
    'UseParallel',true,...
    'Algorithm', 'interior-point',...%'SubproblemAlgorithm', 'factorization', ... 
    ..., %'disp', 'iter',...
    'MaxIterations', 50000);
    [x] = fmincon(fun,x0,[],[],[],[],lb,ub,nonlcon,options);
%     [x] = fminunc(fun,x0);
    c = constraint_test(x,V,Xwb)
%     if c2 > 0.005
%         flag = flag + 1;
%     end
    x_time = reshape(x(1:s,:),[K_stages,2^K_stages])'; %size : 2^K_stages * K_stages
    x_Temp = reshape(x(s+1:end,:),[K_stages,2^K_stages])'; %size : 2^K_stages * K_stages
    expected_Temp = sum(repmat(V', [1 K_stages]).* x_Temp,1);
    expected_time = sum(repmat(V', [1 K_stages]).* x_time,1);
    V_new = association_weights(x,b,V,Xwb);
    V = V_new;
    max_prob = max(V)
    min_prob = min(V)
    x0 = x + pert;
    b=b*alpha;   
    counter = counter - 1;
end
%% Validation of the results


[Vmax,path_opt_index] = max(V);
s = numel(x0)/2;
path_opt = paths(path_opt_index,:);
time = reshape(x(1:s,:),[K_stages,2^K_stages])'; %size : 2^K_stages * K_stages
time_opt = time(path_opt_index,:);
Temp = reshape(x(s+1:end,:),[K_stages,2^K_stages])'; %size : 2^K_stages * K_stages
Temp_opt = Temp(path_opt_index,:);
c = constraint_test(x,V,Xwb);

%% Visualization

t = cell(K_stages,1); mc = cell(K_stages,1); X = Xwb;
for i = 1:K_stages
    t{i} = 0:0.1:time_opt(i);
    mc{i} = zeros(size(t{i}));
    for k = 1:size(t{i},2)
        mc{i}(k) = prediction(Temp_opt(i),path_opt(i),t{i}(k),X);
    end
    X = mc{i}(end);
end

%% Comparison with single stage

T_US0 = []; T_US1 = []; Mwb_f = c;
for T = 30:0.1:70
    if abs(prediction(T,0,sum(time_opt),Xwb) - Mwb_f) < 0.001
        T_US0 = [T_US0 , T];
    end
    if abs(prediction(T,1,sum(time_opt),Xwb) - Mwb_f) < 0.001
        T_US1 = [T_US1 , T];
    end
end
if max(size(T_US1)) == 0
    T_US1 = [30];
end
T_US1 = max(T_US1);
T_US0 = max(T_US0);

tl = 0;
for i = 1:K_stages
    tt = t{i} + tl;
    tl = tt(end);
    figure(1)
    axis square
%     title('Wet-based Moisture Content profile')
%     xlabel('time (min)')
%     ylabel('MC (%)')
    axis ([0 sum(time_opt) 0 100]);
    figure(2)
    axis square
%     title('Temperature profile')
%     xlabel('time (min)')
%     ylabel('Temperature (C)')
    axis ([0 sum(time_opt) 20 80]);
    if rem(path_opt(i),2) == 0
        figure(1)
        f1_1 = plot(tt,100*mc{i},'b');grid on; hold on;
        figure(2)
        f2_1 = plot(tt,Temp_opt(i)*ones(size(t{i})),'b');grid on; hold on;
    else 
        figure(1)
        f1_1 = plot(tt,100*mc{i},'r');grid on; hold on;
        figure(2)
        f2_1 = plot(tt,Temp_opt(i)*ones(size(t{i})),'r');grid on; hold on;
    end
end  
%%

mc_com1 = [];mc_com2 = []; t_comp = 0:0.1:sum(time_opt);
for t_com1 = 0:0.1:sum(time_opt)
    mc_com1 = [mc_com1, prediction(T_US0,0,t_com1,Xwb)];
end
for t_com2 = 0:0.1:sum(time_opt)
    mc_com2 = [mc_com2, prediction(T_US1,1,t_com2,Xwb)];
end
% figure(1)
% f1_3 = plot(t_comp,100*mc_com1,'green');hold on; 
% f1_4 = plot(t_comp,100*mc_com2,'black');
% % f1_2 = plot(t_comp,100*mc_com1,'green');
% % f1_2 = plot(t_comp,100*mc_com2,'black');
% % legend([f1_1,f1_3,f1_4],{'Para-SDM multistage','single stage (Hot air)',...
% %     'single stage (Hot air/Ultrasound)'})
% figure(2)
% f2_3 = plot(t_comp,T_US0*ones(size(t_comp)),'green');hold on; 
% f2_4 = plot(t_comp,T_US1*ones(size(t_comp)),'black');
% % f2_2 = plot(t_comp,T_US0*ones(size(t_comp)),'green');
% % f2_2 = plot(t_comp,T_US1*ones(size(t_comp)),'black');
% % legend([f2_1,f2_3,f2_4] , {'Para-SDM multistage','single stage (Hot air)',...
% %     'single stage (Hot air/Ultrasound)'})
%% Validation (Energy)

p = 1.1;
A = 0.2 * 0.4;
V_air = 1.8;
cp_air = 1.007*1000;
T0 = 20;
Mi = 3.16;
P_us = 168;
gamma = 1/5;
coeff = gamma * p * A * V_air * cp_air;
E_US1 = (P_us + coeff * (T_US1 - T0)) * sum(time_opt);
E_US0 = coeff * (T_US0 - T0) * sum(time_opt);
E_sol = sum(time_opt .* path_opt * P_us + coeff * (Temp_opt - T0) .* time_opt);
