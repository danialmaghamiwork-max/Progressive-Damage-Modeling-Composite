% تابع برای خواندن خواص مواد از فایل اکسل
function materials = read_materials(filename)
    % خواندن فایل اکسل
    data = readtable(filename);
    if height(data) == 0
        error('فایل %s خالی است یا داده‌ای ندارد!', filename);
    end
    materials = struct();
    for i = 1:height(data)
        materials(i).name = string(data.Name{i});
        materials(i).Ex = data.Ex(i) * 1e9;  % تبدیل GPa به Pa
        materials(i).Ey = data.Ey(i) * 1e9;  % تبدیل GPa به Pa
        materials(i).nux = data.nux(i);
        materials(i).Es = data.Es(i) * 1e9;  % تبدیل GPa به Pa
        
        % اضافه کردن استحکام‌ها (MPa به Pa)
        strength_cols = {'Xt', 'Xc', 'Yt', 'Yc', 'S'};
        for j = 1:length(strength_cols)
            col = strength_cols{j};
            if ismember(col, data.Properties.VariableNames) && ~isnan(data.(col)(i))
                materials(i).(col) = data.(col)(i) * 1e6;  % MPa to Pa
            else
                materials(i).(col) = NaN;
            end
        end
        
        % اضافه کردن ضرایب Tsai-Wu با مقیاس جدول (quadratic *1e-18, linear *1e-9)
        tsai_cols_quad = {'Fxx', 'Fyy', 'Fss', 'Fxy'};
        tsai_cols_lin = {'Fx', 'Fy'};
        for j = 1:length(tsai_cols_quad)
            col = tsai_cols_quad{j};
            if ismember(col, data.Properties.VariableNames) && ~isnan(data.(col)(i))
                materials(i).(col) = data.(col)(i) * 1e-18;  % به Pa^{-2}
            else
                materials(i).(col) = NaN;
            end
        end
        for j = 1:length(tsai_cols_lin)
            col = tsai_cols_lin{j};
            if ismember(col, data.Properties.VariableNames) && ~isnan(data.(col)(i))
                materials(i).(col) = data.(col)(i) * 1e-9;  % به Pa^{-1}
            else
                materials(i).(col) = NaN;
            end
        end
    end
end