function Z = TAP(W, a, b, Tlimit)
%TAP 此处显示有关此函数的摘要
%   a and b are row vectors
    mh = rand(1,size(W,2));
    mv = rand(1,size(W,1));
    flag = true;
    mh_pre = rand(1,size(W,2));
    mv_pre = rand(1,size(W,1));
    T = inf;
    T_pre = 0;
    while flag
        abs(T-T_pre)
        if abs(T-T_pre)>Tlimit
            mh_pre = mh;
            mv_pre = mv;
            mh = sigm(b + mv*W - (mv-mv.^2)*W.^2.*(mh-0.5) + ...
                (2*(mh.^2-mh)+1/3).*(((0.5-mv).*(mv-mv.^2))*(W.^3)));
            mv = sigm(a + (W*mh')' - (W.^2*(mh-mh.^2)'.*(mv-0.5)')' + ...
                (2*(mv.^2-mv)+1/3).*(((0.5-mh).*(mh-mh.^2))*(W.^3)'));
        else
            flag = false;
        end
        S_pre = - (mh_pre*log(mh_pre') + (1-mh_pre)*log((1-mh_pre)') + ...
            mv_pre*log(mv_pre') + (1-mv_pre)*log((1-mv_pre)')); 
        T_pre = -S_pre - a*mv_pre' - b*mh_pre' - mv_pre*W*mh_pre' - ...
            0.5*(mv_pre-mv_pre.^2)*W.^2*(mh_pre-mh_pre.^2)' - ...
            2/3*(((mv_pre-mv_pre.^2).*(0.5-mv_pre))*W.^3*((mh_pre-mh_pre.^2).*(0.5-mh_pre))');
        S = - (mh*log(mh') + (1-mh)*log((1-mh)') + mv*log(mv') + (1-mv)*log((1-mv)')); 
        T = -S - a*mv' - b*mh' - mv*W*mh' - 0.5*(mv-mv.^2)*W.^2*(mh-mh.^2)' - ...
            2/3*(((mv-mv.^2).*(0.5-mv))*W.^3*((mh-mh.^2).*(0.5-mh))');
        if isnan(T)
            flag = false;
            T = T_pre;
        end
    end
   
    Z = -T;
end
