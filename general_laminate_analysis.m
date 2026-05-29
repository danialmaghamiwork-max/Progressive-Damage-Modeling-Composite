% general_laminate_analysis.m
% کد اصلی برای تحلیل تنش-کرنش کامپوزیت‌های لایه‌ای عمومی (general laminates)
% این کد بر اساس کدهای قبلی ساخته شده و ماتریس‌های A, B, D را محاسبه می‌کند
% سپس تحت بارهای N و M، کرنش‌های میانی و کرواتورها را حل می‌کند
% و تنش‌ها و کرنش‌ها را برای هر لایه محاسبه می‌کند
% فرض: هسته (core) اختیاری است، اما برای سادگی, بدون هسته در نظر گرفته شده مگر اینکه مشخص شود
% اصلاح: استفاده از متغیرهای invariant U و V برای محاسبه A (برای B و D گسترش خواهد شد)
clear; clc;

% خواندن خواص مواد از فایل اکسل
material_file = 'materials.xlsx';
materials = read_materials(material_file);

% خواندن چیدمان لایه‌ها از فایل اکسل
layup_file = 'laminate_config.xlsx';
[layer_numbers, angles_deg, thickness_mm] = read_laminate_config(layup_file);

% نمایش چیدمان لایه‌ها که خوانده شده
disp('چیدمان لایه‌ها که از فایل خوانده شده:');
for i = 1:length(layer_numbers)
    fprintf('لایه %d: زاویه = %.1f درجه, ضخامت = %.3f mm\n', layer_numbers(i), angles_deg(i), thickness_mm(i));
end

% پرس‌وجو برای هسته (core)
has_core = input('آیا هسته (core) وجود دارد؟ (0 برای خیر, 1 برای بله): ');
if has_core == 1
    core_thickness_mm = input('ضخامت هسته (به میلی‌متر) را وارد کنید: ');
    core_thickness = core_thickness_mm * 1e-3; % تبدیل به متر
else
    core_thickness = 0;
end

% تبدیل ضخامت‌ها به متر
thickness = thickness_mm * 1e-3;

% تعداد لایه‌ها
n_layers = length(angles_deg);

% محاسبه موقعیت‌های z (نسبت به مرکز)
z_bottoms = cumsum([0, thickness(1:end-1)]) - sum(thickness)/2;
z_tops = z_bottoms + thickness;
z_mids = (z_bottoms + z_tops) / 2;

% اگر هسته وجود داشته باشد, موقعیت هسته را تنظیم کنید (فرض وسط)
if has_core
    core_bottom = -core_thickness / 2;
    core_top = core_thickness / 2;
    disp('توجه: برای هسته, موقعیت‌ها تنظیم می‌شود اما کد ساده‌سازی شده است.');
end

% انتخاب مواد برای هر زاویه منحصربه‌فرد
angle_materials = select_materials_for_angles(materials, angles_deg);

% فراخوانی تابع برای محاسبه ماتریس‌های A, B, D با استفاده از U و V
[A, B, D] = calc_abd_matrices(angles_deg, thickness, z_bottoms, z_tops, angle_materials, has_core, core_thickness);

% نمایش ماتریس‌های A, B, D (در GPa*m)
disp('ماتریس A (GPa*m):');
disp(A / 1e9);
disp('ماتریس B (GPa*m^2):');
disp(B / 1e9);
disp('ماتریس D (GPa*m^3):');
disp(D / 1e9);

% ورودی بارهای N و M به صورت جداگانه برای جلوگیری از خطا
disp('ورودی بارهای N (N/m):');
N1 = input('N1: ');
N2 = input('N2: ');
N6 = input('N6: ');
N = [N1; N2; N6];

disp('ورودی بارهای M (Nm/m):');
M1 = input('M1: ');
M2 = input('M2: ');
M6 = input('M6: ');
M = [M1; M2; M6];

loads = [N; M];  % loads 6x1

% محاسبه inv_ABD و strains_curvs
ABD = [A B; B D];
inv_ABD = inv(ABD);
strains_curvs = inv_ABD * loads;
eps0 = strains_curvs(1:3);
kappa = strains_curvs(4:6);

% جمع‌آوری داده برای پلات off-axis
sigma_xx = zeros(n_layers * 2, 1);
sigma_yy = zeros(n_layers * 2, 1);
tau_xy = zeros(n_layers * 2, 1);
z_plot = zeros(n_layers * 2, 1);
idx = 1;

for k = 1:n_layers
    if has_core && abs(z_tops(k) - z_bottoms(k) - core_thickness) < 1e-6
        continue; % پرش هسته
    end
    theta = deg2rad(angles_deg(k));
    selected_material = angle_materials{k};
    [Q_local, ~] = calc_matrices(selected_material.Ex, selected_material.Ey, selected_material.nux, selected_material.Es);
    
    % ماتریس Q_bar
    m = cos(theta); n = sin(theta);
    c2 = m^2; s2 = n^2; c4 = m^4; s4 = n^4; c3s = m^3 * n; cs3 = m * n^3; c2s2 = m^2 * n^2;
    Q11 = Q_local(1,1); Q12 = Q_local(1,2); Q22 = Q_local(2,2); Q66 = Q_local(3,3);
    Qbar11 = Q11*c4 + 2*(Q12 + 2*Q66)*c2s2 + Q22*s4;
    Qbar12 = Q12*(c4 + s4) + (Q11 + Q22 - 4*Q66)*c2s2;
    Qbar16 = (Q11 - Q12 - 2*Q66)*c3s + (Q12 - Q22 + 2*Q66)*cs3;
    Qbar22 = Q11*s4 + 2*(Q12 + 2*Q66)*c2s2 + Q22*c4;
    Qbar26 = (Q11 - Q12 - 2*Q66)*cs3 + (Q12 - Q22 + 2*Q66)*c3s;
    Qbar66 = (Q11 + Q22 - 2*Q12 - 2*Q66)*c2s2 + Q66*(c4 + s4);
    Q_bar = [Qbar11, Qbar12, Qbar16; Qbar12, Qbar22, Qbar26; Qbar16, Qbar26, Qbar66];
    
    % تنش جهانی در bottom و top برای پلات (در Pa)
    eps_global_bottom = eps0 + z_bottoms(k) * kappa;
    eps_global_top = eps0 + z_tops(k) * kappa;
    sigma_global_bottom = Q_bar * eps_global_bottom;
    sigma_global_top = Q_bar * eps_global_top;
    
    % ذخیره برای پلات
    z_plot(idx) = z_bottoms(k);
    sigma_xx(idx) = sigma_global_bottom(1);
    sigma_yy(idx) = sigma_global_bottom(2);
    tau_xy(idx) = sigma_global_bottom(3);
    idx = idx + 1;
    
    z_plot(idx) = z_tops(k);
    sigma_xx(idx) = sigma_global_top(1);
    sigma_yy(idx) = sigma_global_top(2);
    tau_xy(idx) = sigma_global_top(3);
    idx = idx + 1;
end

% پرس‌وجو برای رسم نمودارهای OFF-AXIS
plot_off_axis = input('آیا می‌خواهید نمودارهای OFF-AXIS رسم شوند؟ (1 برای بله, 0 برای خیر): ');
if plot_off_axis == 1
    % رسم نمودارهای جداگانه برای تنش‌های off-axis (در MPa, z در mm)
    figure(1);
    plot(sigma_xx / 1e6, z_plot * 1e3, 'b-', 'LineWidth', 1.5);  % Pa to MPa
    xlabel('Stress \sigma_x (MPa)');
    ylabel('z (mm)');
    title('Off-Axis Stress \sigma_x ');
    grid on;

    figure(2);
    plot(sigma_yy / 1e6, z_plot * 1e3, 'r-', 'LineWidth', 1.5);  % Pa to MPa
    xlabel('Stress \sigma_y (MPa)');
    ylabel('z (mm)');
    title('Off-Axis Stress \sigma_y ');
    grid on;

    figure(3);
    plot(tau_xy / 1e6, z_plot * 1e3, 'g-', 'LineWidth', 1.5);  % Pa to MPa
    xlabel('Stress \tau_{xy} (MPa)');
    ylabel('z (mm)');
    title('Off-Axis Stress \tau_{xy} ');
    grid on;
end

% فراخوانی برای on-axis
calc_onaxis_stress_strain(eps0, kappa, layer_numbers, angles_deg, thickness_mm, z_bottoms, z_tops, angle_materials, has_core, core_thickness);

% پرس‌وجو برای تحلیل progressive damage
do_progressive = input('آیا می‌خواهید تحلیل progressive damage اجرا شود؟ (1 برای بله, 0 برای خیر): ');
if do_progressive == 1
    z_mids = (z_bottoms + z_tops) / 2;
    criterion = input('معیار شکست را انتخاب کنید: 1=Modified Hashin, 2=Hashin, 3=Max Stress, 4=Tsai-Wu: ');
    if criterion == 1
        [FPF, LPF] = progressive_failure_modified_hashin(A, B, D, angles_deg, thickness, z_bottoms, z_tops, z_mids, angle_materials, has_core, core_thickness, layer_numbers);
    elseif criterion == 2
        [FPF, LPF] = progressive_failure_hashin(A, B, D, angles_deg, thickness, z_bottoms, z_tops, z_mids, angle_materials, has_core, core_thickness, layer_numbers);
    elseif criterion == 3
        [FPF, LPF] = progressive_failure_max_stress(A, B, D, angles_deg, thickness, z_bottoms, z_tops, z_mids, angle_materials, has_core, core_thickness, layer_numbers);
    elseif criterion == 4
        [FPF, LPF] = progressive_failure_tsai_wu(A, B, D, angles_deg, thickness, z_bottoms, z_tops, z_mids, angle_materials, has_core, core_thickness, layer_numbers);
    end
end