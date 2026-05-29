function failure_results = hashin_failure(sigma1_layers, sigma2_layers, tau12_layers, Xt_layers, Xc_layers, Yt_layers, Yc_layers, S_layers, layer_numbers)
% Computes the strength ratio R and the failure mode for Hashin criterion
% Inputs: stresses sigma1, sigma2, tau12, strengths Xt, Xc, Yt, Yc, S (all positive)
% Output: min R, and mode string

n_layers = length(sigma1_layers);
failure_results = cell(n_layers, 1);
R_min_layers = inf(n_layers, 1);

disp('تحلیل شکست با معیار Hashin برای هر لایه:');
for k = 1:n_layers
    sigma1 = sigma1_layers(k);
    sigma2 = sigma2_layers(k);
    tau12 = tau12_layers(k);
    Xt = Xt_layers(k);
    Xc = Xc_layers(k);
    Yt = Yt_layers(k);
    Yc = Yc_layers(k);
    S = S_layers(k);  % S = tau_A = S_12 (in-plane shear)
    S_T = Yc / 2;  % S_yz = tau_T = S_23 (transverse shear, approximation Yc / 2)

    if isnan(Xt) || isnan(Xc) || isnan(Yt) || isnan(Yc) || isnan(S)
        failure = false;
        mode = 'No strengths defined';
        R = inf;
    else
        R_candidates = [];
        modes = {};

        if sigma1 >= 0
            a = (sigma1 / Xt)^2 + (tau12 / S)^2;
            if a > 0
                R_ft = 1 / sqrt(a);
                R_candidates = [R_candidates, R_ft];
                modes = [modes, 'Fiber Tension'];
            end
        end

        if sigma1 < 0
            if abs(sigma1) > 0
                R_fc = Xc / abs(sigma1);
                R_candidates = [R_candidates, R_fc];
                modes = [modes, 'Fiber Compression'];
            end
        end

        if sigma2 >= 0
            a = (sigma2 / Yt)^2 + (tau12 / S)^2;
            if a > 0
                R_mt = 1 / sqrt(a);
                R_candidates = [R_candidates, R_mt];
                modes = [modes, 'Matrix Tension'];
            end
        end

        if sigma2 < 0
            a = (sigma2^2 / (4 * S_T^2)) + (tau12^2 / S^2);
            b = ((Yc / (2 * S_T))^2 - 1) * sigma2 / Yc;
            if a > 0
                discriminant = b^2 + 4 * a;
                if discriminant >= 0
                    sqrt_disc = sqrt(discriminant);
                    R1 = (-b + sqrt_disc) / (2 * a);
                    R2 = (-b - sqrt_disc) / (2 * a);
                    R_mc = max([R1, R2]);  % take positive one
                    if R_mc > 0
                        R_candidates = [R_candidates, R_mc];
                        modes = [modes, 'Matrix Compression'];
                    end
                end
            end
        end

        if isempty(R_candidates)
            R = Inf;
            mode = 'No Stress';
        else
            [R, idx] = min(R_candidates);
            mode = modes{idx};
        end
        failure = (R <= 1);
    end

    failure_results{k} = {failure, mode, R};
    R_min_layers(k) = R;

    fprintf('لایه %d: شکست = %s, مد شکست = %s, Strength Ratio (R) = %.4f\n', layer_numbers(k), mat2str(failure), mode, R);
end

if ~isempty(R_min_layers(isfinite(R_min_layers)))
    [min_R_overall, first_failed_idx] = min(R_min_layers);
    fprintf('شکست اولین لایه (FPF) در لایه %d با R = %.4f و مد: %s\n', layer_numbers(first_failed_idx), min_R_overall, failure_results{first_failed_idx}{2});
else
    disp('هیچ لایه‌ای استحکام تعریف‌شده ندارد.');
end
end