function [A, B, D] = calc_abd_matrices(angles_deg, thickness, z_bottoms, z_tops, angle_materials, has_core, core_thickness)
    n_layers = length(angles_deg);
    A = zeros(3,3);
    B = zeros(3,3);
    D = zeros(3,3);
    for k = 1:n_layers
        if has_core && abs(z_tops(k) - z_bottoms(k) - core_thickness) < 1e-6
            continue; % پرش از هسته
        end
        selected_material = angle_materials{k};
        Ex = selected_material.Ex;
        Ey = selected_material.Ey;
        nux = selected_material.nux;
        Es = selected_material.Es;
        [Q_local, ~] = calc_matrices(Ex, Ey, nux, Es);
        U_k = calc_invariantsQ(Q_local); % U برای لایه k
        theta = deg2rad(angles_deg(k));
        cos2 = cos(2*theta);
        cos4 = cos(4*theta);
        sin2 = sin(2*theta);
        sin4 = sin(4*theta);
        dz = thickness(k);
        dz2 = (z_tops(k)^2 - z_bottoms(k)^2) / 2;
        dz3 = (z_tops(k)^3 - z_bottoms(k)^3) / 3;
        % اضافه به A
        A(1,1) = A(1,1) + U_k(1)*dz + U_k(2)*cos2*dz + U_k(3)*cos4*dz;
        A(2,2) = A(2,2) + U_k(1)*dz - U_k(2)*cos2*dz + U_k(3)*cos4*dz;
        A(1,2) = A(1,2) + U_k(4)*dz - U_k(3)*cos4*dz;
        A(2,1) = A(1,2); % symmetric
        A(3,3) = A(3,3) + U_k(5)*dz - U_k(3)*cos4*dz;
        A(1,3) = A(1,3) + (U_k(2)/2)*sin2*dz + U_k(3)*sin4*dz;
        A(3,1) = A(1,3);
        A(2,3) = A(2,3) + (U_k(2)/2)*sin2*dz - U_k(3)*sin4*dz;
        A(3,2) = A(2,3);
        % اضافه به B (مشابه اما با dz2 به جای dz)
        B(1,1) = B(1,1) + U_k(1)*dz2 + U_k(2)*cos2*dz2 + U_k(3)*cos4*dz2;
        B(2,2) = B(2,2) + U_k(1)*dz2 - U_k(2)*cos2*dz2 + U_k(3)*cos4*dz2;
        B(1,2) = B(1,2) + U_k(4)*dz2 - U_k(3)*cos4*dz2;
        B(2,1) = B(1,2);
        B(3,3) = B(3,3) + U_k(5)*dz2 - U_k(3)*cos4*dz2;
        B(1,3) = B(1,3) + (U_k(2)/2)*sin2*dz2 + U_k(3)*sin4*dz2;
        B(3,1) = B(1,3);
        B(2,3) = B(2,3) + (U_k(2)/2)*sin2*dz2 - U_k(3)*sin4*dz2;
        B(3,2) = B(2,3);
        % اضافه به D (مشابه اما با dz3 به جای dz)
        D(1,1) = D(1,1) + U_k(1)*dz3 + U_k(2)*cos2*dz3 + U_k(3)*cos4*dz3;
        D(2,2) = D(2,2) + U_k(1)*dz3 - U_k(2)*cos2*dz3 + U_k(3)*cos4*dz3;
        D(1,2) = D(1,2) + U_k(4)*dz3 - U_k(3)*cos4*dz3;
        D(2,1) = D(1,2);
        D(3,3) = D(3,3) + U_k(5)*dz3 - U_k(3)*cos4*dz3;
        D(1,3) = D(1,3) + (U_k(2)/2)*sin2*dz3 + U_k(3)*sin4*dz3;
        D(3,1) = D(1,3);
        D(2,3) = D(2,3) + (U_k(2)/2)*sin2*dz3 - U_k(3)*sin4*dz3;
        D(3,2) = D(2,3);
    end
end