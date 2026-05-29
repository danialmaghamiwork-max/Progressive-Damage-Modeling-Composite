function calc_onaxis_stress_strain(eps0, kappa, layer_numbers, angles_deg, thickness_mm, z_bottoms, z_tops, angle_materials, has_core, core_thickness)
    % calc_onaxis_stress_strain.m
    % تابع برای محاسبه و نمایش تنش‌ها و کرنش‌های on-axis هر لایه
    % همچنین رسم 6 نمودار: 3 برای تنش on-axis و 3 برای کرنش on-axis بر حسب ضخامت z
    % ورودی‌ها: eps0 (کرنش میانی), kappa (کرواتور), layer_numbers, angles_deg, thickness_mm, z_bottoms, z_tops, angle_materials, has_core, core_thickness

    thickness = thickness_mm * 1e-3; % به متر
    n_layers = length(angles_deg);
    z_mids = (z_bottoms + z_tops) / 2;
    
    % جمع‌آوری داده برای پلات on-axis
    sigma1 = zeros(n_layers * 2, 1); sigma2 = zeros(n_layers * 2, 1); tau12 = zeros(n_layers * 2, 1);
    eps1 = zeros(n_layers * 2, 1); eps2 = zeros(n_layers * 2, 1); gamma12 = zeros(n_layers * 2, 1);
    z_plot = zeros(n_layers * 2, 1);
    idx = 1;
    
    % آرایه‌های جدید برای تنش‌های mid (برای تحلیل شکست)
    sigma1_mid = zeros(n_layers, 1);
    sigma2_mid = zeros(n_layers, 1);
    tau12_mid = zeros(n_layers, 1);
    
    % آرایه‌ها برای استحکام‌ها (لایه‌به‌لایه)
    Xt_layers = NaN(n_layers, 1);
    Xc_layers = NaN(n_layers, 1);
    Yt_layers = NaN(n_layers, 1);
    Yc_layers = NaN(n_layers, 1);
    S_layers = NaN(n_layers, 1);
    
    % آرایه‌ها برای ضرایب Tsai-Wu (اگر از اکسل خوانده شوند)
    Fxx_layers = NaN(n_layers, 1);
    Fyy_layers = NaN(n_layers, 1);
    Fss_layers = NaN(n_layers, 1);
    Fx_layers = NaN(n_layers, 1);
    Fy_layers = NaN(n_layers, 1);
    Fxy_layers = NaN(n_layers, 1);
    
    disp('تنش‌ها و کرنش‌های هر لایه در on-axis:');
    for k = 1:n_layers
        if has_core && abs(z_tops(k) - z_bottoms(k) - core_thickness) < 1e-6
            continue; % پرش از هسته
        end
        theta = deg2rad(angles_deg(k));
        selected_material = angle_materials{k};
        [Q_local, ~] = calc_matrices(selected_material.Ex, selected_material.Ey, selected_material.nux, selected_material.Es);
        
        % کرنش جهانی در bottom, mid, top
        eps_global_bottom = eps0 + z_bottoms(k) * kappa;
        eps_global_mid = eps0 + z_mids(k) * kappa;
        eps_global_top = eps0 + z_tops(k) * kappa;
        
        % ماتریس تبدیل T برای کرنش (engineering notation)
        m = cos(theta); n = sin(theta);
        T_strain = [m^2, n^2, 2*m*n; n^2, m^2, -2*m*n; -m*n, m*n, m^2-n^2];
        
        % کرنش محلی در bottom, mid, top
        eps_local_bottom = T_strain * eps_global_bottom;
        eps_local_mid = T_strain * eps_global_mid;
        eps_local_top = T_strain * eps_global_top;
        
        % تنش محلی در bottom, mid, top
        sigma_local_bottom = Q_local * eps_local_bottom;
        sigma_local_mid = Q_local * eps_local_mid;
        sigma_local_top = Q_local * eps_local_top;
        
        % ذخیره برای پلات (bottom و top)
        z_plot(idx) = z_bottoms(k);
        sigma1(idx) = sigma_local_bottom(1);
        sigma2(idx) = sigma_local_bottom(2);
        tau12(idx) = sigma_local_bottom(3);
        eps1(idx) = eps_local_bottom(1);
        eps2(idx) = eps_local_bottom(2);
        gamma12(idx) = eps_local_bottom(3);
        idx = idx + 1;
        
        z_plot(idx) = z_tops(k);
        sigma1(idx) = sigma_local_top(1);
        sigma2(idx) = sigma_local_top(2);
        tau12(idx) = sigma_local_top(3);
        eps1(idx) = eps_local_top(1);
        eps2(idx) = eps_local_top(2);
        gamma12(idx) = eps_local_top(3);
        idx = idx + 1;
        
        % ذخیره تنش mid برای تحلیل شکست
        sigma1_mid(k) = sigma_local_mid(1);
        sigma2_mid(k) = sigma_local_mid(2);
        tau12_mid(k) = sigma_local_mid(3);
        
        % جمع‌آوری استحکام‌ها و ضرایب
        Xt_layers(k) = selected_material.Xt;
        Xc_layers(k) = selected_material.Xc;
        Yt_layers(k) = selected_material.Yt;
        Yc_layers(k) = selected_material.Yc;
        S_layers(k) = selected_material.S;
        
        Fxx_layers(k) = selected_material.Fxx;
        Fyy_layers(k) = selected_material.Fyy;
        Fss_layers(k) = selected_material.Fss;
        Fx_layers(k) = selected_material.Fx;
        Fy_layers(k) = selected_material.Fy;
        Fxy_layers(k) = selected_material.Fxy;
    end
    
    % پرس‌وجو برای رسم نمودارهای ON-AXIS
    plot_on_axis = input('آیا می‌خواهید نمودارهای ON-AXIS رسم شوند؟ (1 برای بله, 0 برای خیر): ');
    if plot_on_axis == 1
        % رسم نمودارهای on-axis برای تنش‌ها
        figure(4);
        plot(sigma1 / 1e6, z_plot * 1e3, 'b-', 'LineWidth', 1.5);  % Pa to MPa
        xlabel('\sigma_1 (MPa)');
        ylabel('z (mm)');
        title('(on-axis) \sigma_1');
        grid on;
        
        figure(5);
        plot(sigma2 / 1e6, z_plot * 1e3, 'r-', 'LineWidth', 1.5);  % Pa to MPa
        xlabel('\sigma_2 (MPa)');
        ylabel('z (mm)');
        title('(on-axis) \sigma_2');
        grid on;
        
        figure(6);
        plot(tau12 / 1e6, z_plot * 1e3, 'g-', 'LineWidth', 1.5);  % Pa to MPa
        xlabel('\tau_{12} (MPa)');
        ylabel('z (mm)');
        title('(on-axis) \tau_{12}');
        grid on;
        
        % رسم نمودارهای on-axis برای کرنش‌ها
        figure(7);
        plot(eps1, z_plot * 1e3, 'b-', 'LineWidth', 1.5);
        xlabel('\epsilon_1');
        ylabel('z (mm)');
        title('(on-axis) \epsilon_1');
        grid on;
        
        figure(8);
        plot(eps2, z_plot * 1e3, 'r-', 'LineWidth', 1.5);
        xlabel('\epsilon_2');
        ylabel('z (mm)');
        title('(on-axis) \epsilon_2');
        grid on;
        
        figure(9);
        plot(gamma12, z_plot * 1e3, 'g-', 'LineWidth', 1.5);
        xlabel('\gamma_{12}');
        ylabel('z (mm)');
        title('(on-axis) \gamma_{12}');
        grid on;
    end
    
    % سؤال از کاربر برای انتخاب معیارهای شکست (کاربرپسندتر)
    disp('کدام معیارهای شکست را می‌خواهید اجرا کنید؟');
    disp('1: حداکثر تنش محلی');
    disp('2: Tsai-Wu');
    disp('3: Hashin استاندارد');
    disp('4: Hashin اصلاح‌شده');
    disp('می‌توانید چند معیار را انتخاب کنید، مثلاً 1,3 یا همه با 1-4. برای رد، خالی وارد کنید.');
    user_input_str = input('انتخاب شما (مثلاً 1,2,4): ', 's');
    
    if isempty(user_input_str)
        disp('تحلیل شکست رد شد.');
        return;
    end
    
    % تبدیل ورودی به آرایه اعداد
    try
        criteria = str2num(user_input_str); %#ok<ST2NM>
        criteria = unique(criteria(criteria >= 1 & criteria <= 4)); % حذف تکرارها و چک محدوده
        if isempty(criteria)
            disp('ورودی نامعتبر! تحلیل شکست انجام نمی‌شود.');
            return;
        end
    catch
        disp('ورودی نامعتبر! تحلیل شکست انجام نمی‌شود.');
        return;
    end
    
    % اجرای تحلیل شکست برای معیارهای انتخاب‌شده
    if ~isempty(criteria)
        if any(criteria == 1)
            failure_results_max = max_on_axis_stress_failure(sigma1_mid, sigma2_mid, tau12_mid, Xt_layers, Xc_layers, Yt_layers, Yc_layers, S_layers, layer_numbers);
        end
        if any(criteria == 2)
            failure_results_tsai = tsai_wu_failure(sigma1_mid, sigma2_mid, tau12_mid, Xt_layers, Xc_layers, Yt_layers, Yc_layers, S_layers, Fxx_layers, Fyy_layers, Fss_layers, Fx_layers, Fy_layers, Fxy_layers, layer_numbers);
        end
        if any(criteria == 3)
            failure_results_hashin = hashin_failure(sigma1_mid, sigma2_mid, tau12_mid, Xt_layers, Xc_layers, Yt_layers, Yc_layers, S_layers, layer_numbers);
        end
        if any(criteria == 4)
            failure_results_mod_hashin = modified_hashin_failure(sigma1_mid, sigma2_mid, tau12_mid, Xt_layers, Xc_layers, Yt_layers, Yc_layers, S_layers, layer_numbers);
        end
    end
end