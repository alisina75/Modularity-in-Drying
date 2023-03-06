function L = logH(Mwb_final)
    Mwb_desired = 0.07;
    d = (Mwb_final - Mwb_desired)/(Mwb_desired);
    gamma = 100;
    L = 10000 * log((1/gamma) * log(exp(0.001) + exp(gamma * d)));
end