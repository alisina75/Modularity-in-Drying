function  Mwb_out2 = prediction(T2, us2, t2, Mwb_out1)
    % predicting the energy consumption of the stage
%     p = 1.1;
%     A = 0.2 * 0.4;
%     V_air = 1.8;
%     cp_air = 1.007 * 1000;
    %E2 = 200 * t2 * us2 + p * A * V_air * cp_air * T2 * t2;
    % predicting the output moisture content of the stage
    Mi = 3.16;
    T2 = T2 + 273.15;
    if us2 == 1
        k = (-0.005811*T2^2 + 3.9962*T2 - 668.01)/100;
        Meq = (0.1468*T2^2 - 107.27*T2 + 19720)/10000;
        Mdb_out1 = Mwb_out1/(1 - Mwb_out1);
        tin = (1/k) * log ((Mi - Meq)/(Mdb_out1 - Meq));
        t = tin + t2;
        MR2 = exp(-k * t);
        Mdb_out2 = MR2 * (Mi - Meq) + Meq;
        Mwb_out2 = Mdb_out2/(1 + Mdb_out2);
        
    elseif us2 == 0
        k = (0.0074493*T2^2 - 4.5058*T2 + 683.99)/100;
        Meq = (0.2479*T2^2 - 172.09*T2 + 30133)/10000;
        Mdb_out1 = Mwb_out1/(1 - Mwb_out1);
        tin = (1/k) * log ((Mi - Meq)/(Mdb_out1 - Meq));
        t = tin + t2;
        MR2 = exp(-k * t);
        Mdb_out2 = MR2 * (Mi - Meq) + Meq;
        Mwb_out2 = Mdb_out2/(1 + Mdb_out2);
    end
end
        
        