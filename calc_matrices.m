function [Q, S] = calc_matrices(Ex, Ey, nux, Es)
    % محاسبه نسبت پواسون عرضی (nu_y)
    nuy = nux * (Ey / Ex);

    % محاسبه ضریب m
    m = 1 / (1 - nux * nuy);

    % محاسبه ماتریس سختی (Q)
    Qxx = m * Ex;
    Qyy = m * Ey;
    Qxy = m * nuy * Ex;
    Qyx = m * nux * Ey; 
    Qss = Es;

    % ساخت ماتریس سختی (Q)
    Q = [Qxx Qxy 0;
         Qyx Qyy 0;
         0   0   Qss];

    % محاسبه ماتریس نرمی (S)
    Sxx = 1 / Ex;
    Syy = 1 / Ey;
    Sxy = -nuy / Ey;
    Syx = -nux / Ex;
    Sss = 1 / Es;

    % ساخت ماتریس نرمی (S)
    S = [Sxx Sxy 0;
         Syx Syy 0;
         0   0   Sss];
end