function failure_results = max_on_axis_stress_failure(sigma1_layers, sigma2_layers, tau12_layers, Xt_layers, Xc_layers, Yt_layers, Yc_layers, S_layers, layer_numbers)
    % Subroutine for failure analysis using maximum on-axis stress criterion for multiple layers.
    % Includes strength ratio (R) calculation.
    % Inputs:
    %   sigma1_layers: array of sigma1 (on-axis longitudinal stress) for each layer (Pa)
    %   sigma2_layers: array of sigma2 (on-axis transverse stress) for each layer (Pa)
    %   tau12_layers: array of tau12 (on-axis shear stress) for each layer (Pa)
    %   Xt_layers, etc.: arrays of strengths for each layer in Pa (NaN if not defined)
    %   layer_numbers: array of layer numbers for display
    % Outputs:
    %   failure_results: cell array with {failure (logical), mode (string), R_min (double)} for each layer

    n_layers = length(sigma1_layers);
    failure_results = cell(n_layers, 1);
    R_min_layers = inf(n_layers, 1);  % برای ذخیره R_min هر لایه
    
    disp('تحلیل شکست با معیار حداکثر تنش محوری برای هر لایه:');
    for k = 1:n_layers
        sigma_x = sigma1_layers(k);  % sigma1 = sigma_x on-axis
        sigma_y = sigma2_layers(k);  % sigma2 = sigma_y on-axis
        tau_s = tau12_layers(k);     % tau12 = tau_s on-axis
        
        Xt = Xt_layers(k);
        Xc = Xc_layers(k);
        Yt = Yt_layers(k);
        Yc = Yc_layers(k);
        S = S_layers(k);
        
        % اگر استحکام‌ها تعریف نشده، skip
        if isnan(Xt) || isnan(Xc) || isnan(Yt) || isnan(Yc) || isnan(S)
            failure = false;
            mode = 'No strengths defined';
            R_min = inf;
        else
            % محاسبه R برای هر مد
            if sigma_x ~= 0
                if sigma_x > 0
                    R_x = Xt / sigma_x;
                else
                    R_x = Xc / abs(sigma_x);
                end
            else
                R_x = inf;
            end
            
            if sigma_y ~= 0
                if sigma_y > 0
                    R_y = Yt / sigma_y;
                else
                    R_y = Yc / abs(sigma_y);
                end
            else
                R_y = inf;
            end
            
            if tau_s ~= 0
                R_s = S / abs(tau_s);
            else
                R_s = inf;
            end
            
            % پیدا کردن R_min و مد مربوطه
            R_values = [R_x, R_y, R_s];
            modes = {'Fiber failure', 'Matrix failure', 'Shear failure'};
            [R_min, min_idx] = min(R_values);
            failure = (R_min <= 1);
            if failure
                mode = modes{min_idx};
            else
                mode = 'No failure';
            end
        end
        
        failure_results{k} = {failure, mode, R_min};
        R_min_layers(k) = R_min;
        
        % نمایش نتیجه برای لایه فعلی
        fprintf('لایه %d: شکست = %s, مد شکست = %s, Strength Ratio (R_min) = %.4f\n', layer_numbers(k), mat2str(failure), mode, R_min);
    end
    
    % چک کردن شکست کلی (First Ply Failure) بر اساس کوچکترین R_min (فقط لایه‌هایی با R finite)
    valid_R = R_min_layers(isfinite(R_min_layers));
    if ~isempty(valid_R)
        [min_R_overall, first_failed_idx] = min(valid_R);
        fprintf('شکست اولین لایه (FPF) در لایه %d با R = %.4f و مد: %s\n', layer_numbers(first_failed_idx), min_R_overall, failure_results{first_failed_idx}{2});
    else
        disp('هیچ لایه‌ای استحکام تعریف‌شده ندارد.');
    end
end