function [c, ceq] = nonlin_const1(x0,V,Xwb)
    s = numel(x0)/2;
    K_stages = log(numel(V))/log(2);
    time0 = reshape(x0(1:s,:),[K_stages,2^K_stages])'; %size : 2^K_stages * K_stages
    Temp0 = reshape(x0(s+1:end,:),[K_stages,2^K_stages])'; %size : 2^K_stages * K_stages
    K_stages = log(numel(V))/log(2);
    num_paths = numel(V);
    paths = dec2bin(0:2^K_stages-1) - '0'; %size : 2^K_stages * K_stages 
    %% final constraints
    Mwb_final = zeros(num_paths,1);
    for i=1:num_paths 
        if imag(path_prediction(paths(i,:),Temp0(i,:),time0(i,:),Xwb)) ~= 0
            flag1 = 1;
        end
        Mwb_final(i,1) = path_prediction(paths(i,:),Temp0(i,:),time0(i,:),Xwb);
    end
    %c = total_time - time_lim;
    c = [];
    %c = [c1;c2];
    ceq = Mwb_final - 0.075 * ones(size(Mwb_final));
    %ceq = [];
end