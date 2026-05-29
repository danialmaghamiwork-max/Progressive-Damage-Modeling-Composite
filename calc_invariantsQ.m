function U = calc_invariantsQ(Q)
    % تابع برای محاسبه متغیرهای ناوردا U
    % ورودی: ماتریس سختی Q
    % خروجی: U = [U1Q; U2Q; U3Q; U4Q; U5Q]
    Qxx = Q(1,1); Qyy = Q(2,2); Qxy = Q(1,2); Qss = Q(3,3);
    U1Q = (1/8) * (3 * Qxx + 3 * Qyy + 2 * Qxy + 4 * Qss);
    U2Q = (1/2) * (Qxx - Qyy);
    U3Q = (1/8) * (Qxx + Qyy - 2 * Qxy - 4 * Qss);
    U4Q = (1/8) * (Qxx + Qyy + 6 * Qxy - 4 * Qss);
    U5Q = (1/8) * (Qxx + Qyy - 2 * Qxy + 4 * Qss);
    U = [U1Q; U2Q; U3Q; U4Q; U5Q];
end