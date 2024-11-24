clear all;
clc;

Luminance_range = logspace(log10(1),log10(1000),20);
refresh_rate_range = [24,36,48,60,90,120];

% Assumed VRR Display contrast - luminance relation
b = log10(0.008);
k = (b - log10(0.0025)) / log10(100);
contrast_function = @(RR, Lum) ((144 - RR) / (144 - 40)) * 10 ^ ((-k * log10(Lum) + b));
RR_inverse_function = @(contrast, Lum) 144 - contrast ./ (10 .^ ((-k .* log10(Lum) + b))) .* (144 - 40);

csf_elaTCSF_model = CSF_elaTCSF();

calculate_need = 0;
plot_need = 1;

T_frequency_range = linspace(0,20,50);
Width = 62.7;
Height = 37.8;
Peak_sensitivity_list = [];

if calculate_need == 1
    for Luminance_index = 1:length(Luminance_range)
        Luminance_value = Luminance_range(Luminance_index);
        sensitivity_list = [];
        for t_frequency_index = 1:length(T_frequency_range)
            t_frequency_value = T_frequency_range(t_frequency_index);
            csf_pars = struct('s_frequency', 0, 't_frequency', t_frequency_value, 'orientation', 0, ...
                'luminance', Luminance_value, 'width', Width, 'height', Height, 'eccentricity', 0);
            S = csf_elaTCSF_model.sensitivity_rect(csf_pars);
            sensitivity_list(end+1) = S;
        end
        peak_sensitivity = max(sensitivity_list);
        Peak_sensitivity_list(end+1) = peak_sensitivity;
        writematrix(Peak_sensitivity_list, 'application_vrr_range_peak_sensitivity_list.csv');
    end
else
    Peak_sensitivity_list = readmatrix('application_vrr_range_peak_sensitivity_list.csv');
end

ha = tight_subplot(1, 2, [.13 .09],[.15 .02],[.1 .03]);
set(gcf, 'Position', [100, 100, 700, 320]);
axes(ha(1));
Contrast_Ticks = [0.001,0.01];
Contrast_Threshold_Curve = 1 ./ Peak_sensitivity_list;

x_fill = [flip(Luminance_range), Luminance_range];
y_fill = [ones(size(Contrast_Threshold_Curve)).*min(Contrast_Ticks), Contrast_Threshold_Curve];
fill(x_fill, y_fill, [1 0.8 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
hold on;

hh = [];
for RR_index = 1:length(refresh_rate_range)
    RR_value = refresh_rate_range(RR_index);
    contrast_list = [];
    for Lum_index = 1:length(Luminance_range)
        Luminance_value = Luminance_range(Lum_index);
        contrast = contrast_function(RR_value, Luminance_value);
        contrast_list(end+1) = contrast;
    end
    hh(end+1) = plot(Luminance_range, contrast_list, 'LineWidth', 2, 'DisplayName', [num2str(RR_value) ' Hz']);
    hold on;
end
plot(Luminance_range, Contrast_Threshold_Curve, 'r', 'LineWidth', 2, 'LineStyle', '--');
x_center = 10^mean([log10(min(Luminance_range)), log10(max(Luminance_range))]);
text(x_center, min(Contrast_Threshold_Curve), sprintf('Contrast Threshold\npredicted by elaTCSF'), ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontSize', 10);

xlabel('Luminance (cd/m^2)', 'FontSize', 12);
ylabel('\textbf{\textnormal{VRR Contrast }} $C = \Delta L / \overline{L} $', 'FontSize', 12, 'Interpreter', 'latex');
xlim([min(Luminance_range), max(Luminance_range)]);
ylim([min(Contrast_Ticks), max(Contrast_Ticks)]);
set(gca, 'XScale', 'log', 'XTick', [0.1, 1, 10, 100, 1000], 'XTickLabel', [0.1, 1, 10, 100, 1000]);
set(gca, 'YScale', 'log', 'YTick', Contrast_Ticks, 'YTickLabel', Contrast_Ticks);
legend (hh, 'Location', 'best','NumColumns', 2);

axes(ha(2));
RR_ticks = [24, 36, 48, 60, 90, 120, 144];
RR_list = RR_inverse_function(Contrast_Threshold_Curve, Luminance_range);
x_fill = [Luminance_range, flip(Luminance_range)];
y_fill = [RR_list, ones(size(RR_list)).*144];
fill(x_fill, y_fill, [0.8 1 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
hold on;
plot(Luminance_range, RR_list, 'g', 'LineWidth', 2);
text(100, 100, sprintf('Flicker\nInvisible'), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', 12);
set(gca, 'XScale', 'log', 'XTick', [0.1, 1, 10, 100, 1000], 'XTickLabel', [0.1, 1, 10, 100, 1000]);
xlabel('Luminance (cd/m^2)', 'FontSize', 12);
ylabel('Refresh Rate (Hz)', 'FontSize', 12);
xlim([min(Luminance_range), max(Luminance_range)]);
ylim([min(RR_ticks), max(RR_ticks)]);
yticks(RR_ticks);
yticklabels(RR_ticks);
