log_coeff =...
[-18.601440, -18.751220, -18.968840, -19.133310;...
 -18.565540, -18.715990, -18.934000, -19.099460;...
 -18.530450, -18.680260, -18.899830, -19.066350;...
 -18.485350, -18.634280, -18.855230, -19.022510]

logTe = [1.301230, 1.477320, 1.699170, 1.845300];
logNe = [18.301030, 18.698970, 19.000000, 19.301030];

f_sub = zeros(4,4);

f_sub(1:2, 1:2) = log_coeff(2:3, 2:3)

% for i = 1:2
% for j = 1:2
	
% end
% end