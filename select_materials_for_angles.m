% تابع برای انتخاب مواد برای هر زاویه منحصربه‌فرد
function [angle_materials] = select_materials_for_angles(materials, angles_deg)
    % پیدا کردن زوایای منحصربه‌فرد
    unique_angles = unique(angles_deg);
    angle_materials = cell(length(angles_deg), 1); % آرایه سلولی برای مواد هر لایه
    
    % انتخاب ماده برای هر زاویه منحصربه‌فرد
    for i = 1:length(unique_angles)
        current_angle = unique_angles(i);
        disp(['انتخاب ماده برای زاویه ', num2str(current_angle), ' درجه:']);
        layers_with_angle = find(angles_deg == current_angle);
        
        % نمایش لیست مواد
        disp('لیست مواد موجود:');
        for j = 1:length(materials)
            fprintf('%d. %s\n', j, materials(j).name);
        end
        
        % دریافت شماره ماده از کاربر
        prompt = sprintf('لطفاً شماره ماده را برای زاویه %d درجه وارد کنید (1 تا %d): ', current_angle, length(materials));
        selected_index = input(prompt);
        if selected_index < 1 || selected_index > length(materials)
            error('شماره ماده نامعتبر است! لطفاً شماره بین 1 تا %d وارد کنید.', length(materials));
        end
        
        % اختصاص ماده به همه لایه‌های با زاویه فعلی
        selected_material = materials(selected_index);
        for j = 1:length(layers_with_angle)
            angle_materials{layers_with_angle(j)} = selected_material;
        end
    end
end

