function failure_results = modified_hashin_failure(sigma1_layers, sigma2_layers, tau12_layers, Xt_layers, Xc_layers, Yt_layers, Yc_layers, S_layers, layer_numbers)
% Computes the strength ratio R and the failure mode for modified Hashin criterion
% Inputs: stresses sigma1, sigma2, tau12, strengths Xt, Xc, Yt, Yc, S (all positive)
% Output: min R, and mode string

n_layers = length(sigma1_layers);
failure_results = cell(n_layers, 1);
R_min_layers = inf(n_layers, 1);

disp('تحلیل شکست با معیار Modified Hashin برای هر لایه:');
for k = 1:n_layers
    sigma1 = sigma1_layers(k);
    sigma2 = sigma2_layers(k);
    tau12 = tau12_layers(k);
    Xt = Xt_layers(k);
    Xc = Xc_layers(k);
    Yt = Yt_layers(k);
    Yc = Yc_layers(k);
    S = S_layers(k);

    if isnan(Xt) || isnan(Xc) || isnan(Yt) || isnan(Yc) || isnan(S)
        failure = false;
        mode = 'No strengths defined';
        R = inf;
    else
        R_candidates = [];
        modes = {};

        % Fiber Breakage in Tension (sigma1 > 0)
        if sigma1 > 0
            a = (sigma1 / Xt)^2 + (tau12 / S)^2;
            if a > 0
                R_fbt = 1 / sqrt(a);
                R_candidates = [R_candidates, R_fbt];
                modes = [modes, 'Fiber Breakage Tension'];
            end
        end

        % Fiber Buckling in Compression (sigma1 < 0)
        if sigma1 < 0
            if abs(sigma1) > 0
                R_fbc = Xc / abs(sigma1);
                R_candidates = [R_candidates, R_fbc];
                modes = [modes, 'Fiber Buckling Compression'];
            end
        end

        % Matrix Tension (sigma2 > 0)
        if sigma2 > 0
            a = (sigma2 / Yt)^2 + (tau12 / S)^2;
            if a > 0
                R_mt = 1 / sqrt(a);
                R_candidates = [R_candidates, R_mt];
                modes = [modes, 'Matrix Tension'];
            end
        end

        % Matrix Compression (sigma2 < 0)
        if sigma2 < 0
            a = (sigma2 / Yc)^2 + (tau12 / S)^2;
            if a > 0
                R_mc = 1 / sqrt(a);
                R_candidates = [R_candidates, R_mc];
                modes = [modes, 'Matrix Compression'];
            end
        end

        % Fiber-Matrix Shearing in Compression (sigma1 < 0)
        if sigma1 < 0
            a = (abs(sigma1) / Xc)^2 + (tau12 / S)^2;
            if a > 0
                R_fms = 1 / sqrt(a);
                R_candidates = [R_candidates, R_fms];
                modes = [modes, 'Fiber-Matrix Shearing Compression'];
            end
        end

        if isempty(R_candidates)
            R = Inf;
            mode = 'No Stress';
        else
            [R, idx] = min(R_candidates);
            mode = modes{idx};
        end
        failure = (R < 1); % تغییر به <1 برای سازگاری با >1
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