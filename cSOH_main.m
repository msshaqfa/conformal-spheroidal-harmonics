% Fast ellipsoidal conformal map (FECM) for genus-0 closed surfaces
% and spheroidal harmonic (SOH) decomposition.
% This code combines two different works as detailed below.
% 
% If you use this code in your work, please cite the following papers:
%
% [1] G. P. T. Choi, 
%     "Fast ellipsoidal conformal and quasi-conformal parameterization of genus-0 closed surfaces".
%     Journal of Computational and Applied Mathematics, 447, 115888, 2024.
% 
% [2] M. Shaqfa, W. M. van Rees,
%   "Spheroidal harmonics for generalizing the morphological decomposition of closed parametric
%   surfaces". Construction and Building Materials, 2024.
% 
% Copyright (c) 2024, Mahmoud Shaqfa and Gary P. T. Choi

clear; clc; close all;
addpath('code');
addpath('stlTools/');
addpath('Input_geometry/');
addpath('Input_geometry/data');
%% Input surface


% filename = "buddha.mat"; % works n_max = 70
% filename = "vaselion.mat"; % works n_max = 70
% filename = "bimba.mat"; % works for low resolution n_max = 35
% filename = "max_ref_poisson.stl";
% filename = "bulldog.mat"; % works n_max = 30
% filename = "nefertiti.mat"; % works n_max = 30

[~, baseName, file_extension] = fileparts(filename);
if file_extension == ".mat"
    load(filename);
elseif file_extension == ".stl"
    [v, f, ~, ~] = stlRead(filename);
end

mapping_type = "conformal";

max_n = 20;
rec_max_n = max_n;

rec_refinement = 6; % max ~6, otherwise it gets very large.

plot_mesh(v,f);
view([-70 10]);
title('Input surface');

%% Surface registration (Alignment and placement) & spheroidal sizing

[new_v, centroid, Rots] = register_surface(v, false);
[aa, cc, spheroid_type, R_y] = fit_spheroids(v, false);
foci =  sqrt(abs(aa^2 - cc^2));
new_v = (R_y * new_v')';

if spheroid_type == "oblate" % for oblates
    zeta = acosh(aa / foci);
else % for prolates and spheres
    zeta = asinh(aa / foci);
end

scatter3(new_v(:, 1), new_v(:, 2), new_v(:, 3), 'kx', 'DisplayName', 'New boundary')
box on; axis equal;
title("Boundary fitting for surface registratrion")
xlabel('x')
ylabel('y')
zlabel('z')

hold on

func = @(x_,y_,z_) (x_.^2 + y_.^2)/ (aa^2) + z_.^2 / (cc^2) - 1;
fimplicit3( func,'EdgeColor', 'none', 'FaceAlpha', .5)
hold off

%% Surface Parameterization (Conformal or Quasi-conformal mapping)
if mapping_type == "conformal"
    map = ellipsoidal_conformal_map(new_v, f, aa, aa, cc);
end

plot_mesh(map, f);
view([-70 10]);
title('Ellipsoidal conformal parameterization');

angle_distortion = angle_distortion(new_v, f, map);
area_distortion = area_distortion(new_v, f, map);

figure;
histogram(angle_distortion,-180:1:180);
xlim([-180 180])
title('Angle Distortion');
xlabel('Angle difference (degree)')
ylabel('Number of angles')
set(gca,'FontSize',12);

figure;
histogram(area_distortion,-5:0.1:5);
xlim([-5 5])
title('Area Distortion');
xlabel('log(final area/initial area)')
ylabel('Number of faces')
set(gca,'FontSize',12);

%% Basis functions and decomposition

% Spheroidal coordinates
[thetas, phis] = cart2spheroid(map, foci, zeta, spheroid_type);

% Spheroidal basis functions
D_mat = spheroidal_harmonic_basis(max_n, thetas, phis, spheroid_type);

% Expand the harmonics (find Fourier weights)
qm_k = D_mat\new_v;

%% Shape descriptors and the fractal dimension
% Shape descriptors (2-norm) for frequency accumulates at a certain
% frequency degree..
Dl = zeros([3, max_n+1]);
for k_ = 1:3
    for n_ = 1:max_n
        for m_ = -n_:1:n_
            Dl(k_, n_) = Dl(k_, n_) + (real(qm_k(n_^2 + n_ + m_ + 1, k_)))^2 ...
                + (imag(qm_k(n_^2 + n_ + m_ + 1, k_)))^2;
        end
        Dl(k_, n_) = sqrt(Dl(k_, n_));
    end
    Dl(k_, :) = Dl(k_, :)/Dl(k_, 1);
end

figure
subplot(2,3,1)
x_temp = 1:length(Dl(1, :));
loglog(x_temp(2:end), Dl(1, 2:end), 'LineWidth', 2)
title('The shape descriptors')
xlabel('Freqency index (k)')
ylabel('Normalized amplitude (Dx)')
grid on

subplot(2,3,2)
loglog(x_temp(2:end), Dl(2, 2:end), 'LineWidth', 2)
title('The shape descriptors')
xlabel('Freqency index (k)')
ylabel('Normalized amplitude (Dy)')
grid on

subplot(2,3,3)
loglog(x_temp(2:end), Dl(3, 2:end), 'LineWidth', 2)
title('The shape descriptors')
xlabel('Freqency index (k)')
ylabel('Normalized amplitude (Dz)')
grid on

subplot(2,3,4:6)
loglog(x_temp(2:end), sqrt(Dl(1, 2:end).^2 + ...
    Dl(2, 2:end).^2 + Dl(3, 2:end).^2), 'LineWidth', 2)
title('The shape descriptors')
xlabel('Freqency index (k)')
ylabel('Normalized amplitude (Dr)')


%% Surface reconstruction
[rec_v, rec_f] = icosahedron_sphere(rec_refinement, false);

% Scale it to a spheroid (ignore the distortion)
rec_v(:, 1) = rec_v(:, 1) * aa;
rec_v(:, 2) = rec_v(:, 2) * aa;
rec_v(:, 3) = rec_v(:, 3) * cc;

[rec_thetas, rec_phis] = cart2spheroid(rec_v, foci, zeta, spheroid_type);


rec_D_mat = spheroidal_harmonic_basis(rec_max_n, rec_thetas, ...
    rec_phis, spheroid_type);

rec_v = real(rec_D_mat * qm_k(1:(rec_max_n+1)^2, :));

rec_v_old = rec_v;

% To restore the surface registration
rec_v = (R_y' * rec_v')' * Rots';

rec_v = rec_v + centroid;

paraview_patch(rec_v, rec_f)
fname_out = strcat('./output/rec_n_', num2str(rec_max_n), '_', mapping_type, '_',baseName, '.stl');
if endsWith(fname_out, '.mat', 'IgnoreCase', true)
    fname_out = strrep(fname_out, '.mat', '.stl');
end
stlWrite(fname_out, rec_f, rec_v);