load('acd96_c.mat')

log_coeff = reshape(log_coeff_json, [24, 30, 6]);

log_coeff(12, 15, 1) = -4;

Te_sample = linspace(min(log_temperature), max(log_temperature), 100);
Ne_sample = linspace(min(log_density), max(log_density), 100);
Coeff = zeros(100,100);

k = 1;

for Te_it = 1:length(Te_sample)
% for Ne_it = 1:length(Ne_sample)
Te = Te_sample(Te_it);
% Ne = Ne_sample(Ne_it)
Coeff(:, Te_it) = interp2(log_temperature, log_density, squeeze(log_coeff(:, :, k)), Te, Ne_sample, 'spline');
% end
end

% size(Coeff)

s = surf(Te_sample, Ne_sample, Coeff, 'FaceAlpha', 0.5, 'EdgeColor', 'none');

hold on

surf(log_temperature, log_density, squeeze(log_coeff(:, :, k)), 'FaceAlpha', 0);

% k_test = [4];
% Te_test = 50;
Ne_test = 1e19;


% for k_it = 1:length(k_test)
% for Te_it = 1:length(Te_test)
% for Ne_it = 1:length(Ne_test)

% k = k_test(k_it);
% Te = Te_test(Te_it);
% Ne = Ne_test(Ne_it);

	% surf(log_temperature, log_density, squeeze(log_coeff(:, :, k)));

% 	Vq = interp2(log_temperature, log_density, squeeze(log_coeff(:, :, k)), log10(Te), log10(Ne), 'spline');

% 	fprintf('log_coeff(%d, %f, %e) = %f \n', k, Te, Ne, Vq)

% end
% end
% end

clear Ne Te k Vq Te_sample Ne_sample