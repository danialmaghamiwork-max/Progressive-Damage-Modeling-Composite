function [FPF_load, LPF_load] = progressive_failure_modified_hashin(A_init, B_init, D_init, angles_deg, thickness, z_bottoms, z_tops, z_mids, angle_materials, has_core, core_thickness, layer_numbers)
    % progressive_failure_modified_hashin.m
    % محاسبه FPF و LPF با progressive damage modeling با استفاده از Modified Hashin
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
    for k = 1:n_layers
        Xt_layers(k) = materials_current{k}.Xt;
        Xc_layers(k) = materials_current{k}.Xc;
        Yt_layers(k) = materials_current{k}.Yt;
        Yc_layers(k) = materials_current{k}.Yc;
        S_layers(k) = materials_current{k}.S;
    end
    
    % آرایه برای ذخیره ترتیب شکست لایه‌ها
    failed_order = [];
    
    disp('شروع تحلیل progressive damage با Modified Hashin...');
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
            if has_core && abs(thickness(k) - core_thickness) < 1e-6 || failed_layers(k)
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
        
        % چک شکست با Modified Hashin
        failure_results = modified_hashin_failure(sigma1_mid, sigma2_mid, tau12_mid, Xt_layers, Xc_layers, Yt_layers, Yc_layers, S_layers, layer_numbers);
        
        % پیدا کردن لایه‌های جدید شکست‌خورده
        new_failures = false;
        for k = 1:n_layers
            if ~failed_layers(k) && failure_results{k}{1}  % failure = true
                mode = failure_results{k}{2};
                fprintf('لایه %d با مد شکست %s شناسایی شد.\n', layer_numbers(k), mode);
                mat_k = materials_current{k};
                if strcmp(mode, 'Fiber Breakage Tension')
                    mat_k.Ex = mat_k.Ex / 10000;
                    mat_k.Ey = mat_k.Ey / 10000;
                    mat_k.nux = mat_k.nux / 10000;
                    mat_k.Es = mat_k.Es / 10000;
                    % به‌روزرسانی دستی ماتریس Q برای حذف اثرات پواسون
                    nuy = mat_k.nux * (mat_k.Ey / mat_k.Ex);
                    m_initial = 1 / (1 - mat_k.nux * nuy);
                    Q_initial = [mat_k.Ex * m_initial, mat_k.Ex * nuy * m_initial, 0; ...
                                 mat_k.Ey * mat_k.nux * m_initial, mat_k.Ey * m_initial, 0; ...
                                 0, 0, mat_k.Es];
                    mat_k.Q_local = [100000 / Q_initial(1,1), 100000 / Q_initial(1,2), 0; ...
                                     100000 / Q_initial(2,1), 100000 / Q_initial(2,2), 0; ...
                                     0, 0, 100000 / Q_initial(3,3)];
                elseif strcmp(mode, 'Fiber Buckling Compression')
                    mat_k.Ex = mat_k.Ex / 10000;
                    mat_k.Ey = mat_k.Ey / 10000;
                    mat_k.nux = mat_k.nux / 10000;
                    mat_k.Es = mat_k.Es / 10000;
                    % به‌روزرسانی دستی ماتریس Q برای حذف اثرات پواسون
                    nuy = mat_k.nux * (mat_k.Ey / mat_k.Ex);
                    m_initial = 1 / (1 - mat_k.nux * nuy);
                    Q_initial = [mat_k.Ex * m_initial, mat_k.Ex * nuy * m_initial, 0; ...
                                 mat_k.Ey * mat_k.nux * m_initial, mat_k.Ey * m_initial, 0; ...
                                 0, 0, mat_k.Es];
                    mat_k.Q_local = [100000 / Q_initial(1,1), 100000 / Q_initial(1,2), 0; ...
                                     100000 / Q_initial(2,1), 100000 / Q_initial(2,2), 0; ...
                                     0, 0, 100000 / Q_initial(3,3)];
                elseif strcmp(mode, 'Matrix Tension')
                    mat_k.Ey = mat_k.Ey / 10000;
                    mat_k.nux = mat_k.nux / 10000;
                    % E_x و G_s بدون تغییر باقی می‌مونن
                    % به‌روزرسانی دستی ماتریس Q
                    nuy = mat_k.nux * (mat_k.Ey / mat_k.Ex);
                    m_initial = 1 / (1 - mat_k.nux * nuy);
                    Q_initial = [mat_k.Ex * m_initial, mat_k.Ex * nuy * m_initial, 0; ...
                                 mat_k.Ey * mat_k.nux * m_initial, mat_k.Ey * m_initial, 0; ...
                                 0, 0, mat_k.Es];
                    mat_k.Q_local = [mat_k.Ex, 100000 / Q_initial(2,2), 100000 / Q_initial(1,2); ...
                                     100000 / Q_initial(2,1), 100000 / Q_initial(2,2), 0; ...
                                     0, 0, mat_k.Es];
                elseif strcmp(mode, 'Matrix Compression')
                    mat_k.Ey = mat_k.Ey / 10000;
                    mat_k.nux = mat_k.nux / 10000;
                    % E_x و G_s بدون تغییر باقی می‌مونن
                    % به‌روزرسانی دستی ماتریس Q
                    nuy = mat_k.nux * (mat_k.Ey / mat_k.Ex);
                    m_initial = 1 / (1 - mat_k.nux * nuy);
                    Q_initial = [mat_k.Ex * m_initial, mat_k.Ex * nuy * m_initial, 0; ...
                                 mat_k.Ey * mat_k.nux * m_initial, mat_k.Ey * m_initial, 0; ...
                                 0, 0, mat_k.Es];
                    mat_k.Q_local = [mat_k.Ex, 100000 / Q_initial(2,2), 100000 / Q_initial(1,2); ...
                                     100000 / Q_initial(2,1), 100000 / Q_initial(2,2), 0; ...
                                     0, 0, mat_k.Es];
                elseif strcmp(mode, 'Fiber-Matrix Shearing Compression')
                    mat_k.nux = mat_k.nux / 10000;
                    mat_k.Es = mat_k.Es / 10000;
                    % E_x و E_y بدون تغییر باقی می‌مونن
                    nuy = mat_k.nux * (mat_k.Ey / mat_k.Ex);
                    m_initial = 1 / (1 - mat_k.nux * nuy);
                    Q_initial = [mat_k.Ex * m_initial, mat_k.Ex * nuy * m_initial, 0; ...
                                 mat_k.Ey * mat_k.nux * m_initial, mat_k.Ey * m_initial, 0; ...
                                 0, 0, mat_k.Es];
                    mat_k.Q_local = [mat_k.Ex, 100000 / Q_initial(1,2), 0; ...
                                     100000 / Q_initial(2,1), mat_k.Ey, 0; ...
                                     0, 0, 100000 / mat_k.Es];
                end
                materials_current{k} = mat_k;
                failed_layers(k) = true;
                new_failures = true;
                failed_order = [failed_order, layer_numbers(k)];  % اضافه کردن لایه به ترتیب شکست
                fprintf('لایه %d شکست در lambda=%.2f با مد: %s\n', layer_numbers(k), lambda, mode);
                if isnan(FPF_load)
                    FPF_load = lambda;
                end
                
                % حلقه داخلی برای چک شکست‌های جدید تو همون lambda
                inner_new_failures = true;
                while inner_new_failures
                    inner_new_failures = false;
                    % به‌روزرسانی ABD
                    [A, B, D] = calc_abd_matrices(angles_deg, thickness, z_bottoms, z_tops, materials_current, has_core, core_thickness);
                    ABD = [A, B; B, D];
                    if det(A) < 1e-10
                        break;
                    end
                    inv_ABD = inv(ABD);
                    eps_kappa = inv_ABD * current_load;
                    eps0 = eps_kappa(1:3);
                    kappa = eps_kappa(4:6);
                    
                    % محاسبه دوباره تنش‌ها
                    sigma1_mid = zeros(n_layers, 1);
                    sigma2_mid = zeros(n_layers, 1);
                    tau12_mid = zeros(n_layers, 1);
                    for inner_k = 1:n_layers
                        if has_core && abs(thickness(inner_k) - core_thickness) < 1e-6 || failed_layers(inner_k)
                            continue;
                        end
                        theta = deg2rad(angles_deg(inner_k));
                        mat_inner = materials_current{inner_k};
                        if isfield(mat_inner, 'Q_local')
                            Q_local = mat_inner.Q_local; % استفاده از Q_local دستی اگه تعریف شده
                        else
                            [Q_local, ~] = calc_matrices(mat_inner.Ex, mat_inner.Ey, mat_inner.nux, mat_inner.Es);
                        end
                        eps_global_mid = eps0 + z_mids(inner_k) * kappa;
                        T_strain = [cos(theta)^2, sin(theta)^2, 2*cos(theta)*sin(theta); ...
                                    sin(theta)^2, cos(theta)^2, -2*cos(theta)*sin(theta); ...
                                    -cos(theta)*sin(theta), cos(theta)*sin(theta), cos(theta)^2 - sin(theta)^2];
                        eps_local_mid = T_strain * eps_global_mid;
                        sigma_local_mid = Q_local * eps_local_mid;
                        sigma1_mid(inner_k) = sigma_local_mid(1);
                        sigma2_mid(inner_k) = sigma_local_mid(2);
                        tau12_mid(inner_k) = sigma_local_mid(3);
                    end
                    
                    % چک دوباره شکست
                    failure_results_new = modified_hashin_failure(sigma1_mid, sigma2_mid, tau12_mid, Xt_layers, Xc_layers, Yt_layers, Yc_layers, S_layers, layer_numbers);
                    
                    for inner_k = 1:n_layers
                        if ~failed_layers(inner_k) && failure_results_new{inner_k}{1}
                            mode = failure_results_new{inner_k}{2};
                            fprintf('شکست زنجیره‌ای: لایه %d با مد %s در lambda=%.2f\n', layer_numbers(inner_k), mode, lambda);
                            mat_inner = materials_current{inner_k};
                            if strcmp(mode, 'Fiber Breakage Tension')
                                mat_inner.Ex = mat_inner.Ex / 10000;
                                mat_inner.Ey = mat_inner.Ey / 10000;
                                mat_inner.nux = mat_inner.nux / 10000;
                                mat_inner.Es = mat_inner.Es / 10000;
                                nuy = mat_inner.nux * (mat_inner.Ey / mat_inner.Ex);
                                m_initial = 1 / (1 - mat_inner.nux * nuy);
                                Q_initial = [mat_inner.Ex * m_initial, mat_inner.Ex * nuy * m_initial, 0; ...
                                             mat_inner.Ey * mat_inner.nux * m_initial, mat_inner.Ey * m_initial, 0; ...
                                             0, 0, mat_inner.Es];
                                mat_inner.Q_local = [100000 / Q_initial(1,1), 100000 / Q_initial(1,2), 0; ...
                                                     100000 / Q_initial(2,1), 100000 / Q_initial(2,2), 0; ...
                                                     0, 0, 100000 / Q_initial(3,3)];
                            elseif strcmp(mode, 'Fiber Buckling Compression')
                                mat_inner.Ex = mat_inner.Ex / 10000;
                                mat_inner.Ey = mat_inner.Ey / 10000;
                                mat_inner.nux = mat_inner.nux / 10000;
                                mat_inner.Es = mat_inner.Es / 10000;
                                nuy = mat_inner.nux * (mat_inner.Ey / mat_inner.Ex);
                                m_initial = 1 / (1 - mat_inner.nux * nuy);
                                Q_initial = [mat_inner.Ex * m_initial, mat_inner.Ex * nuy * m_initial, 0; ...
                                             mat_inner.Ey * mat_k.nux * m_initial, mat_inner.Ey * m_initial, 0; ...
                                             0, 0, mat_inner.Es];
                                mat_inner.Q_local = [100000 / Q_initial(1,1), 100000 / Q_initial(1,2), 0; ...
                                                     100000 / Q_initial(2,1), 100000 / Q_initial(2,2), 0; ...
                                                     0, 0, 100000 / Q_initial(3,3)];
                            elseif strcmp(mode, 'Matrix Tension')
                                mat_inner.Ey = mat_inner.Ey / 10000;
                                mat_inner.nux = mat_inner.nux / 10000;
                                % E_x و G_s بدون تغییر باقی می‌مونن
                                nuy = mat_inner.nux * (mat_inner.Ey / mat_inner.Ex);
                                m_initial = 1 / (1 - mat_inner.nux * nuy);
                                Q_initial = [mat_inner.Ex * m_initial, mat_inner.Ex * nuy * m_initial, 0; ...
                                             mat_inner.Ey * mat_inner.nux * m_initial, mat_inner.Ey * m_initial, 0; ...
                                             0, 0, mat_inner.Es];
                                mat_inner.Q_local = [mat_inner.Ex, 100000 / Q_initial(2,2), 100000 / Q_initial(1,2); ...
                                                     100000 / Q_initial(2,1), 100000 / Q_initial(2,2), 0; ...
                                                     0, 0, mat_inner.Es];
                            elseif strcmp(mode, 'Matrix Compression')
                                mat_inner.Ey = mat_inner.Ey / 10000;
                                mat_inner.nux = mat_inner.nux / 10000;
                                % E_x و G_s بدون تغییر باقی می‌مونن
                                nuy = mat_inner.nux * (mat_inner.Ey / mat_inner.Ex);
                                m_initial = 1 / (1 - mat_inner.nux * nuy);
                                Q_initial = [mat_inner.Ex * m_initial, mat_inner.Ex * nuy * m_initial, 0; ...
                                             mat_inner.Ey * mat_inner.nux * m_initial, mat_inner.Ey * m_initial, 0; ...
                                             0, 0, mat_inner.Es];
                                mat_inner.Q_local = [mat_inner.Ex, 100000 / Q_initial(2,2), 100000 / Q_initial(1,2); ...
                                                     100000 / Q_initial(2,1), 100000 / Q_initial(2,2), 0; ...
                                                     0, 0, mat_inner.Es];
                            elseif strcmp(mode, 'Fiber-Matrix Shearing Compression')
                                mat_inner.nux = mat_inner.nux / 10000;
                                mat_inner.Es = mat_inner.Es / 10000;
                                % E_x و E_y بدون تغییر باقی می‌مونن
                                nuy = mat_inner.nux * (mat_inner.Ey / mat_inner.Ex);
                                m_initial = 1 / (1 - mat_inner.nux * nuy);
                                Q_initial = [mat_inner.Ex * m_initial, mat_inner.Ex * nuy * m_initial, 0; ...
                                             mat_inner.Ey * mat_inner.nux * m_initial, mat_inner.Ey * m_initial, 0; ...
                                             0, 0, mat_inner.Es];
                                mat_inner.Q_local = [mat_inner.Ex, 100000 / Q_initial(1,2), 0; ...
                                                     100000 / Q_initial(2,1), mat_inner.Ey, 0; ...
                                                     0, 0, 100000 / mat_inner.Es];
                            end
                            materials_current{inner_k} = mat_inner;
                            failed_layers(inner_k) = true;
                            failed_order = [failed_order, layer_numbers(inner_k)];
                            inner_new_failures = true;
                        end
                    end
                end
            end
        end
        
        % اگر شکست جدید، گام برگرد و با گام کوچکتر امتحان کن
        if new_failures
            lambda = lambda - lambda_step;
            lambda_step = lambda_step / 10;  % کوچکتر کن
            if lambda_step < 1e-3  % حداقل دقت
                lambda_step = 1e-3;
            end
            continue;
        end
        
        % چک اگر همه شکست
        if all(failed_layers | (abs(thickness - core_thickness) < 1e-6))  % هسته رو حساب نکن
            LPF_load = lambda;
            break;
        end
    end
    
    disp(['FPF Load (N/m): ', num2str(FPF_load)]);
    disp(['LPF Load (N/m): ', num2str(LPF_load)]);
    disp(['ترتیب شکست لایه‌ها: ', num2str(failed_order)]);
end