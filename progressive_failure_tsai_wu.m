function [FPF_load, LPF_load] = progressive_failure_tsai_wu(A_init, B_init, D_init, angles_deg, thickness, z_bottoms, z_tops, z_mids, angle_materials, has_core, core_thickness, layer_numbers)
    % progressive_failure_tsai_wu.m
    % محاسبه FPF و LPF با progressive damage modeling با استفاده از Tsai-Wu
    % ورودی‌ها: ماتریس‌های اولیه A B D، زوایا، ضخامت‌ها، موقعیت z، مواد، هسته، شماره لایه‌ها
    % خروجی: FPF_load و LPF_load (در N/m برای N_x)
    
    n_layers = length(angles_deg);
    % کپی از مواد برای تغییر (degrade)
    materials_current = angle_materials;
    failed_layers = false(n_layers, 1);
    
    % بار base: uniaxial tension N_x = 1 N/m, بقیه صفر
    load_base = [1; 0; 0; 0; 0; 0];  % [N_x, N_y, N_xy, M_x, M_y, M_xy]
    
    % گام اولیه lambda (مقیاس بار، می‌تونی تنظیم کنی)
    lambda_step = 1e3;  % شروع با گام بزرگ، می‌تونی کوچک کن برای دقت
    lambda = 0;
    FPF_load = NaN;
    LPF_load = NaN;
    
    % آرایه برای استحکام‌ها (از مواد استخراج)
    Xt_layers = zeros(n_layers, 1);
    Xc_layers = zeros(n_layers, 1);
    Yt_layers = zeros(n_layers, 1);
    Yc_layers = zeros(n_layers, 1);
    S_layers = zeros(n_layers, 1);
    Fxx_layers = zeros(n_layers, 1);
    Fyy_layers = zeros(n_layers, 1);
    Fss_layers = zeros(n_layers, 1);
    Fx_layers = zeros(n_layers, 1);
    Fy_layers = zeros(n_layers, 1);
    Fxy_layers = zeros(n_layers, 1);
    for k = 1:n_layers
        Xt_layers(k) = materials_current{k}.Xt;
        Xc_layers(k) = materials_current{k}.Xc;
        Yt_layers(k) = materials_current{k}.Yt;
        Yc_layers(k) = materials_current{k}.Yc;
        S_layers(k) = materials_current{k}.S;
        Fxx_layers(k) = materials_current{k}.Fxx;
        Fyy_layers(k) = materials_current{k}.Fyy;
        Fss_layers(k) = materials_current{k}.Fss;
        Fx_layers(k) = materials_current{k}.Fx;
        Fy_layers(k) = materials_current{k}.Fy;
        Fxy_layers(k) = materials_current{k}.Fxy;
    end
    
    % آرایه برای ذخیره ترتیب شکست لایه‌ها
    failed_order = [];
    
    disp('شروع تحلیل progressive damage با Tsai-Wu...');
    while true
        lambda = lambda + lambda_step;
        current_load = lambda * load_base;
        
        % بروزرسانی ABD با مواد فعلی
        [A, B, D] = calc_abd_matrices(angles_deg, thickness, z_bottoms, z_tops, materials_current, has_core, core_thickness);
        ABD = [A, B; B, D];
        if det(A) < 1e-10  % ظرفیت extensional صفر
            LPF_load = lambda - lambda_step;
            break;
        end
        inv_ABD = inv(ABD);
        eps_kappa = inv_ABD * current_load;
        eps0 = eps_kappa(1:3);
        kappa = eps_kappa(4:6);
        
        % محاسبه تنش mid برای چک شکست
        sigma1_mid = zeros(n_layers, 1);
        sigma2_mid = zeros(n_layers, 1);
        tau12_mid = zeros(n_layers, 1);
        for k = 1:n_layers
            if has_core && abs(thickness(k) - core_thickness) < 1e-6
                continue;  % پرش هسته
            end
            if failed_layers(k)
                continue;
            end
            theta = deg2rad(angles_deg(k));
            mat_k = materials_current{k};
            [Q_local, ~] = calc_matrices(mat_k.Ex, mat_k.Ey, mat_k.nux, mat_k.Es);
            eps_global_mid = eps0 + z_mids(k) * kappa;
            T_strain = [cos(theta)^2, sin(theta)^2, 2*cos(theta)*sin(theta); ...
                        sin(theta)^2, cos(theta)^2, -2*cos(theta)*sin(theta); ...
                        -cos(theta)*sin(theta), cos(theta)*sin(theta), cos(theta)^2 - sin(theta)^2];
            eps_local_mid = T_strain * eps_global_mid;
            sigma_local_mid = Q_local * eps_local_mid;
            sigma1_mid(k) = sigma_local_mid(1);
            sigma2_mid(k) = sigma_local_mid(2);
            tau12_mid(k) = sigma_local_mid(3);
        end
        
        % چک شکست با Tsai-Wu
        failure_results = tsai_wu_failure(sigma1_mid, sigma2_mid, tau12_mid, Xt_layers, Xc_layers, Yt_layers, Yc_layers, S_layers, Fxx_layers, Fyy_layers, Fss_layers, Fx_layers, Fy_layers, Fxy_layers, layer_numbers);
        
        % پیدا کردن لایه‌های جدید شکست‌خورده
        new_failures = false;
        for k = 1:n_layers
            if ~failed_layers(k) && failure_results{k}{1}  % failure = true
                mode = failure_results{k}{2};
                fprintf('لایه %d شکست در lambda=%.2f با مد: %s\n', layer_numbers(k), lambda, mode);
                mat_k = materials_current{k};
                % برای Tsai-Wu, dominant stress
                sigma1 = sigma1_mid(k);
                sigma2 = sigma2_mid(k);
                tau12 = tau12_mid(k);
                if abs(sigma1) > abs(sigma2) && abs(sigma1) > abs(tau12)
                    mat_k.Ex = mat_k.Ex / 10000;
                    mat_k.nux = mat_k.nux / 10000;
                    mat_k.Es = mat_k.Es / 10000;
                elseif abs(sigma2) > abs(sigma1) && abs(sigma2) > abs(tau12)
                    mat_k.Ey = mat_k.Ey / 10000;
                    mat_k.nux = mat_k.nux / 10000;
                    mat_k.Es = mat_k.Es / 10000;
                else
                    mat_k.Es = mat_k.Es / 10000;
                end
                materials_current{k} = mat_k;
                failed_layers(k) = true;
                new_failures = true;
                failed_order = [failed_order, layer_numbers(k)];
                if isnan(FPF_load)
                    FPF_load = lambda;
                end
            end
        end
        
        if new_failures
            lambda = lambda - lambda_step;
            lambda_step = lambda_step / 10;
            if lambda_step < 1e-3
                lambda_step = 1e-3;
            end
            continue;
        end
        
        if all(failed_layers | (abs(thickness - core_thickness) < 1e-6))
            LPF_load = lambda;
            break;
        end
    end
    
    disp(['FPF Load (N/m): ', num2str(FPF_load)]);
    disp(['LPF Load (N/m): ', num2str(LPF_load)]);
    disp(['ترتیب شکست لایه‌ها: ', num2str(failed_order)]);
end