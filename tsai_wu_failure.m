function failure_results = tsai_wu_failure(sigma1_layers, sigma2_layers, tau12_layers, Xt_layers, Xc_layers, Yt_layers, Yc_layers, S_layers, Fxx_layers, Fyy_layers, Fss_layers, Fx_layers, Fy_layers, Fxy_layers, layer_numbers)
    % Subroutine for failure analysis using Tsai-Wu criterion for multiple layers.
    % Includes failure index and strength ratio (R) calculation.
    % Inputs: similar to max_stress, plus F coefficients arrays (NaN if not defined in excel)
    % Outputs: cell array with {failure (logical), mode (string), R (double), failure_index (double)} for each layer

    n_layers = length(sigma1_layers);
    failure_results = cell(n_layers, 1);
    R_layers = inf(n_layers, 1);  % برای ذخیره R هر لایه
    FI_layers = zeros(n_layers, 1);  % برای ذخیره شاخص شکست
    
    disp('تحلیل شکست با معیار Tsai-Wu برای هر لایه:');
    for k = 1:n_layers
        sigma_x = sigma1_layers(k);  % sigma1 = sigma_x on-axis (Pa)
        sigma_y = sigma2_layers(k);  % sigma2 = sigma_y on-axis (Pa)
        tau_s = tau12_layers(k);     % tau12 = tau_s on-axis (Pa)
        
        % چک اگر ضرایب از اکسل تعریف شده باشند
        use_predefined_F = ~isnan(Fxx_layers(k)) && ~isnan(Fyy_layers(k)) && ~isnan(Fss_layers(k)) && ...
                           ~isnan(Fx_layers(k)) && ~isnan(Fy_layers(k)) && ~isnan(Fxy_layers(k));
        
        if use_predefined_F
            Fxx = Fxx_layers(k);
            Fyy = Fyy_layers(k);
            Fss = Fss_layers(k);
            Fx = Fx_layers(k);
            Fy = Fy_layers(k);
            Fxy = Fxy_layers(k);
        else
            % محاسبه از استحکام‌ها اگر تعریف نشده
            Xt = Xt_layers(k);
            Xc = Xc_layers(k);
            Yt = Yt_layers(k);
            Yc = Yc_layers(k);
            S = S_layers(k);
            
            if isnan(Xt) || isnan(Xc) || isnan(Yt) || isnan(Yc) || isnan(S)
                failure = false;
                mode = 'No strengths defined';
                R = inf;
                FI = 0;
            else
                Fxx = 1 / (Xt * Xc);
                Fyy = 1 / (Yt * Yc);
                Fss = 1 / (S ^ 2);
                Fx = 1/Xt - 1/Xc;
                Fy = 1/Yt - 1/Yc;
                Fxy = -0.5 * sqrt(Fxx * Fyy);  % تقریب استاندارد
            end
        end
        
        % محاسبه شاخص شکست (FI)
        FI = Fxx * sigma_x^2 + Fyy * sigma_y^2 + Fss * tau_s^2 + ...
             2 * Fxy * sigma_x * sigma_y + Fx * sigma_x + Fy * sigma_y;
        FI_layers(k) = FI;
        
        % محاسبه نسبت استحکام R
        A = Fxx * sigma_x^2 + 2 * Fxy * sigma_x * sigma_y + Fyy * sigma_y^2 + Fss * tau_s^2;
        B = Fx * sigma_x + Fy * sigma_y;
        
        if A == 0
            if B == 0
                R = inf;  % بدون بار
            elseif B < 0
                R = -1 / B;
            else
                R = inf;
            end
        else
            discriminant = B^2 + 4 * A;
            if discriminant < 0
                R = inf;  % بدون شکست واقعی
            else
                sqrt_disc = sqrt(discriminant);
                R1 = (-B + sqrt_disc) / (2 * A);
                R2 = (-B - sqrt_disc) / (2 * A);
                % انتخاب ریشه مثبت
                if R1 * R2 < 0
                    R = max(R1, R2);
                elseif R1 > 0 && R2 > 0
                    R = min(R1, R2);
                else
                    R = inf;
                end
            end
        end
        R_layers(k) = R;
        
        failure = (FI >= 1);
        if failure
            mode = 'Tsai-Wu interactive failure';
        else
            mode = 'No failure';
        end
        
        failure_results{k} = {failure, mode, R, FI};
        
        % نمایش نتیجه برای لایه فعلی
        fprintf('لایه %d: شکست = %s, مد شکست = %s, Strength Ratio (R) = %.4f, Failure Index (FI) = %.4f\n', ...
                layer_numbers(k), mat2str(failure), mode, R, FI);
    end
    
    % چک کردن شکست کلی (First Ply Failure) بر اساس کوچکترین R (فقط لایه‌هایی با R finite)
    valid_R = R_layers(isfinite(R_layers));
    if ~isempty(valid_R)
        [min_R_overall, first_failed_idx] = min(valid_R);
        fprintf('شکست اولین لایه (FPF) در لایه %d با R = %.4f و مد: %s\n', layer_numbers(first_failed_idx), min_R_overall, failure_results{first_failed_idx}{2});
    else
        disp('هیچ لایه‌ای استحکام یا ضریب تعریف‌شده ندارد.');
    end
end