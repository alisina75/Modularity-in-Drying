function V_new = association_weights(x0,b,V,Xwb)
    s = numel(x0)/2;
    K_stages = log(numel(V))/log(2);
    time0 = reshape(x0(1:s,:),[K_stages,2^K_stages])'; %size : 2^K_stages * K_stages
    Temp0 = reshape(x0(s+1:end,:),[K_stages,2^K_stages])'; %size : 2^K_stages * K_stages
    K_stages = log(numel(V))/log(2);
    num_paths = numel(V);
    paths = dec2bin(0:2^K_stages-1) - '0'; %size : 2^K_stages * K_stages 
    %% Energy consumption 
    p = 1.1;
    A = 0.2 * 0.4;
    V_air = 1.8;
    cp_air = 1.007*1000;
    T0 = 20;
    P_us = 168;
    gamma = 1/5;
    coeff = gamma * p * A * V_air * cp_air;
    Temp_dif = Temp0 - T0 * ones(size(Temp0));
    E_HA = coeff * Temp_dif .* time0;
    E_US = P_us * paths .* time0;
    E = sum(E_US + E_HA,2);
    E = E/1000;
    %% final constraints
    Mwb_final = zeros(num_paths,1);
    total_time = sum(time0,2);
    for i=1:num_paths 
        Mwb_final(i,1) = path_prediction(paths(i,:),Temp0(i,:),time0(i,:),Xwb);
    end
    Mwb_p = MC_penalty(Mwb_final);
    D_act = E + Mwb_p;
    D = D_act - min(D_act);
    V_new = (exp(-b*D)./sum(exp(-b*D)))';
end