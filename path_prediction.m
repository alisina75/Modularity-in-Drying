function  Mwb_f = path_prediction(path_row,T_row,t_row,Xwb)
    K_stages = numel(path_row);
    Mwb_in = Xwb;
    for j = 1:K_stages
        Mwb_out = prediction(T_row(1,j),path_row(1,j),t_row(1,j),Mwb_in);
        Mwb_in = Mwb_out;
    end
    Mwb_f = Mwb_in;       
end
