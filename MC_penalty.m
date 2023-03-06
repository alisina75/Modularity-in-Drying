function L = MC_penalty(Mwb_final)
    Mwb_desired = 0.085;
    d1 = (Mwb_final - Mwb_desired)/(Mwb_desired);
    gamma = 100;
    L1 = 1000 * (log((1/gamma) * log(exp(0.001) + exp(gamma * d1)))...
        - log((1/gamma) * log(exp(0.001) + exp(gamma * -1))));
%     d2 = 0.14 - Mwb_final;
%     L2 = 100 * (log((1/gamma) * log(exp(0.001) + exp(gamma * d2)))...
%         - log((1/gamma) * log(exp(0.001) + exp(gamma * -1))));
%     L = L1 + L2;
    L = L1;
end