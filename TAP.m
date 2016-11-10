function Z = TAP(W, a, b)
%TAP 此处显示有关此函数的摘要
%   a and b are row vectors
    mh = zeros(1,size(W,2));
    mv = zeros(1,size(W,1));
    flag = true;
    mh_pre = rand(1,size(W,2));
    mv_pre = rand(1,size(W,1));
    i = 0;
    min_mh = inf;
    min_mv = inf;
    while flag
        if abs(norm(mh-mh_pre))>1.2*min_mh || abs(norm(mv-mv_pre))>1.2*min_mv || i <50000 
            mh_pre = mh;
            mv_pre = mv;
            mh = sigm(b + mv*W - (mv-mv.^2)*W.^2.*(mh-0.5));
            mv = sigm(a + (W*mh')' - (W.^2*(mh-mh.^2)'.*(mv-0.5)')');
        else
            flag = false;
        end
        if i > 30000
            if abs(norm(mh-mh_pre))< min_mh 
                min_mh = abs(norm(mh-mh_pre));
            end
            if abs(norm(mv-mv_pre))< min_mv 
                min_mv = abs(norm(mv-mv_pre));
            end
        end
        i = i+1;
    end
    S = - (mh*log(mh') + (1-mh)*log((1-mh)') + mv*log(mv') + (1-mv)*log((1-mv)')); 
    T = -S - a*mv' - b*mh' - mv*W*mh' - 0.5*(mv-mv.^2)*W.^2*(mh-mh.^2)';
    Z = -T;
end
