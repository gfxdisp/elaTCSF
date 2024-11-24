% elaTCSF
clear all;
clc;
% addpath('../../models/');
% addpath('../../utils/');
% csf_model = CSF_elaTCSF_16();
csf_elaTCSF_model = CSF_elaTCSF();
csf_model = csf_elaTCSF_model;

% IEC filter
% Define the numerator and denominator coefficients of the transfer function
num = [0.041661, 44.758, 2715.6, 29839, 0];
den = [1, 196.32, 11781, 534820, 3505380];

% Define the transfer function
H = tf(num, den);

% Define the frequency range for the plot
f = logspace(-1, 2.1, 1000); % Frequency range from 10^-2 to 10^4 Hz
% f = logspace(-2, 4, 1000); % Frequency range from 10^-2 to 10^4 Hz

omega = 2 * pi * f; % Convert frequency to angular frequency
s = 1./omega;

% Calculate the frequency response
[mag, phase] = bode(H, omega);
mag = squeeze(mag);
mag_norm = mag./max(mag);


% Bodington
f_bod = @(f) 0.01254 - 0.0007571.*f - 0.00004007.*(f.^2) + 0.000006757.*(f.^3) -...
    0.00000023306.*(f.^4) + 0.000000002958.*(f.^5);
M_bod = (1./(f_bod(f)));
M_bod_norm = M_bod./max(M_bod);

% Plot the amplitude response
% figure,
clf;
gap = [0.06 ,0.04];
marg_h = [0.17, 0.33]; % [lower upper]
marg_w = [0.07, 0.02]; % [left right]
% [ha, axh] = tight_subplot(1, 2, gap, marg_h, marg_w);
ha = tight_subplot(1, 2, gap, marg_h, marg_w);
% set(gcf, 'Position', [100, 100, 1500, 800]);

set(gcf, 'defaultAxesFontName', 'Arial');
set(gcf, 'defaultTextFontName', 'Arial');
% set(gcf, 'defaultTextFontSize', 9);
FontSize = 12;

fig = gcf;
fig_size = [5 2.3]*2;
fig.Units = 'inches';
fig.Position = [2, 2, fig_size];
fig.PaperUnits = 'inches';
fig.PaperSize = fig_size;

%% sub plot (luminance)
axes(ha(1));
loglog(f, mag_norm, '--k', 'LineWidth', 1); % Grey line
labels{1} = 'IEC Light-flickermeter';

% title('Amplitude response of chosen filters');

hold on,
loglog(f, M_bod_norm, '--r', 'LineWidth', 1); % Grey line
labels{2} = 'Bodington et al. 2016';

Luminance_range = logspace(log10(1000),log10(1),4);
Area = 1000;
l_cols = lines(numel(Luminance_range));


for ll = 1:length(Luminance_range)
    csf_pars = struct('s_frequency', 0, 't_frequency', f, 'orientation', 0, ...
        'luminance', Luminance_range(ll), 'area', Area, 'eccentricity', 0);
    S = csf_model.sensitivity(csf_pars);
    if ll == 1
        max_factor = max(S);
    end
    S_norm = S./max_factor;
    loglog(f, S_norm, 'Color', l_cols(ll, :), 'LineWidth', 2);
    labels{ll+2} = sprintf('elaTCSF (%g cd/m^2)',  Luminance_range(ll));
end


% xlabel({'Temp. freq. [Hz]', '(a)'});
xlabel({'Temporal Frequency [Hz]', '(a)'});
ylabel('Relative flicker sensitivity');

grid on;
% xlim([10^-2 10^4]);

xticks(10.^(-2:1:4));
% xticklabels([]);

ylim([10^-3.3 10^0]);
yticks(10.^(-3:1:0));


legend(labels, 'Position', [0.1325    0.6897    0.4411    0.2795],...
    'Box','off', 'FontSize', FontSize );

% Adjust plot properties to match the figure style
% set(gca, 'FontSize', 12);
% set(gca, 'LineWidth', 1.5);


%% sub plot (eccentricity)
axes(ha(2));
loglog(f, mag_norm, '--k', 'LineWidth', 1); % Grey line
labels{1} = 'IEC Light-flickermeter';

% title('Amplitude response of chosen filters');

hold on,
loglog(f, M_bod_norm, '--r', 'LineWidth', 1); % Grey line
labels{2} = 'Bodington et al. 2016';

eccentricity_range = [0, 10, 20, 30];
e_cols = lines(length(eccentricity_range));
Area = 100;

for ee = 1:length(eccentricity_range)
    csf_pars = struct('s_frequency', 0, 't_frequency', f, 'orientation', 0, ...
        'luminance',1000, 'area', Area, 'eccentricity', eccentricity_range(ee));
    S = csf_model.sensitivity(csf_pars);
    if ee == 1
        max_factor = max(S);
    end
    S_norm = S./max_factor;
    loglog(f, S_norm, 'Color', e_cols(ee, :), 'LineWidth', 2);
    labels{ee+2} = sprintf('elaTCSF (%g^o)',  eccentricity_range(ee));
end


grid on;
% xlabel({'Temp. freq. [Hz]', '(b)'});
xlabel({'Temporal Frequency [Hz]', '(a)'});
% ylabel('Relative flicker sensitivity');
xticks(10.^[-2:1:4]);

ylim([10^-3.3 10^0]);
yticks(10.^(-3:1:0));
yticklabels([]);


legend(labels, 'Position', [0.5554    0.6931    0.4260    0.2795],...
    'Box','off', 'FontSize', FontSize);


%% figure printing

% Find all text objects in the figure
textObjects = findall(gcf, 'Type', 'text');
set(textObjects, 'FontSize', FontSize);

axesObjects = findall(gcf, 'Type', 'axes');
set(axesObjects, 'FontSize', FontSize);

if 0
    name = 'plot_elatcf_flicker_app';
     print(fig, '-vector', name,'-dpdf');
    print(fig, '-vector', name,'-dpng', '-r600');
end
