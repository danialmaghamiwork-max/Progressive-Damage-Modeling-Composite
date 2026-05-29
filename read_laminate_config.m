% تابع برای خواندن چیدمان لایه‌ها از فایل اکسل با پشتیبانی از تقارن و تکرار
function [layer_numbers, angles_deg, thickness_mm] = read_laminate_config(filename)
    % تابع برای خواندن چیدمان لایه‌ها از فایل اکسل با پشتیبانی از تقارن و تکرار
    % ورودی: نام فایل اکسل
    % خروجی: شماره لایه‌ها (خودکار)، زوایا (درجه)، ضخامت‌ها (میلی‌متر)
    % فرض: فایل اکسل شامل ستون‌های Angle (degrees)، NumLayers، Thickness (mm)، و Symmetry (s یا ns یا T یا nT یا خالی)
    
    % خواندن فایل اکسل
    data = readtable(filename, 'VariableNamingRule', 'preserve');
    
    % چک کردن ستون‌های ضروری
    required_columns = {'Angle (degrees)', 'NumLayers', 'Thickness (mm)'};
    if ~all(ismember(required_columns, data.Properties.VariableNames))
        error('فایل اکسل باید ستون‌های Angle (degrees)، NumLayers، و Thickness (mm) داشته باشد!');
    end
    
    % استخراج داده‌ها
    angles_deg_input = data.("Angle (degrees)");
    num_layers_input = fix(data.NumLayers); % اطمینان از عدد صحیح
    thickness_mm_input = data.("Thickness (mm)");
    if any(~isnumeric(angles_deg_input)) || any(~isnumeric(num_layers_input)) || any(num_layers_input <= 0) || ...
       any(~isnumeric(thickness_mm_input)) || any(thickness_mm_input <= 0)
        error('زوایا، تعداد لایه‌ها، و ضخامت‌ها باید عددی و مثبت باشند!');
    end
    
    % استخراج تقارن/تکرار (اختیاری)
    if ismember('Symmetry', data.Properties.VariableNames)
        symmetry_input = data.Symmetry;
    else
        symmetry_input = repmat("", size(angles_deg_input)); % پیش‌فرض بدون تقارن
    end
    
    % یافتن مقدار Symmetry معتبر (فرض: فقط یکی غیرخالی، یا همه یکسان)
    non_empty_sym = symmetry_input(~cellfun(@isempty, symmetry_input));
    if ~isempty(non_empty_sym)
        unique_sym = unique(non_empty_sym);
        if length(unique_sym) > 1
            error('فقط یک مقدار Symmetry معتبر مجاز است (همه ردیف‌ها باید خالی یا یکسان باشند)!');
        end
        sym = strtrim(lower(unique_sym{1})); % تبدیل به حروف کوچک
    else
        sym = ""; % بدون تقارن
    end
    
    % گسترش زوایا و ضخامت‌ها بر اساس تعداد لایه‌ها (چیدمان پایه)
    base_angles_deg = [];
    base_thickness_mm = [];
    for i = 1:length(angles_deg_input)
        base_angles_deg = [base_angles_deg, repmat(angles_deg_input(i), 1, num_layers_input(i))];
        base_thickness_mm = [base_thickness_mm, repmat(thickness_mm_input(i), 1, num_layers_input(i))];
    end
    
    % اعمال تکرار و تقارن بر اساس Symmetry برای کل چیدمان
    angles_deg = base_angles_deg;
    thickness_mm = base_thickness_mm;
    
    if ~isempty(sym)
        if strcmp(sym, 's') % تقارن ساده: کل چیدمان + برعکس آن
            angles_deg = [angles_deg, fliplr(angles_deg)];
            thickness_mm = [thickness_mm, fliplr(thickness_mm)];
        elseif endsWith(sym, 's') && length(sym) > 1 % ns: تکرار (n-1) بار + تقارن
            repeat_count_str = sym(1:end-1);
            repeat_count = str2double(repeat_count_str);
            if isnan(repeat_count) || repeat_count < 1 || mod(repeat_count, 1) ~= 0
                error('فرمت ns اشتباه است! n باید عدد صحیح مثبت باشد.');
            end
            for r = 1:repeat_count - 1
                angles_deg = [angles_deg, base_angles_deg];
                thickness_mm = [thickness_mm, base_thickness_mm];
            end
            angles_deg = [angles_deg, fliplr(angles_deg)];
            thickness_mm = [thickness_mm, fliplr(thickness_mm)];
        elseif strcmp(sym, 't') % T: بدون تغییر (فقط خواندن)
            % هیچ کاری انجام نشود
        elseif endsWith(sym, 't') && length(sym) > 1 % nT: تکرار (n-1) بار بدون تقارن
            repeat_count_str = sym(1:end-1);
            repeat_count = str2double(repeat_count_str);
            if isnan(repeat_count) || repeat_count < 1 || mod(repeat_count, 1) ~= 0
                error('فرمت nT اشتباه است! n باید عدد صحیح مثبت باشد.');
            end
            for r = 1:repeat_count - 1
                angles_deg = [angles_deg, base_angles_deg];
                thickness_mm = [thickness_mm, base_thickness_mm];
            end
        else
            error('مقدار Symmetry نامعتبر است! باید s, ns, T, nT یا خالی باشد.');
        end
    end
    
    % تولید شماره لایه‌ها به‌صورت خودکار
    layer_numbers = 1:length(angles_deg);
end