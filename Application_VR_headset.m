clear all;
clc;

csf_elaTCSF_model = CSF_elaTCSF();

calculate_need = 0;
plot_need = 1;

contrast_list = [0.1]; % This value is influenced by the scaling factor (as mentioned in the text, varying across datasets). Please conduct further measurements to determine its precise value before use.
peak_luminance_visio_pro = 5000;
peak_luminance_holoLens_2 = 500;
peak_luminance_quest_3 = 100;

if calculate_need == 1
    % Assume The persistence = 0.1
    % Luminance_average = 0.1 * Luminance_Peak
    width = 100;
    height = 100;
    luminance_list_visionpro = logspace(log10(10),log10(10000),20);
    VR_headset_RR_matrix_visionpro = zeros(length(contrast_list), length(luminance_list_visionpro));
    for contrast_index = 1:length(contrast_list)
        contrast = contrast_list(contrast_index);
        find_sensitivity = 1 / contrast;
        for luminance_index = 1:length(luminance_list_visionpro)
            luminance = luminance_list_visionpro(luminance_index);
            bs_func = @(omega) - csf_elaTCSF_model.sensitivity_rect(struct('s_frequency', 0, 't_frequency', omega, 'orientation', 0, ...
                'luminance', luminance.*0.1, 'width', width, 'height', height, 'eccentricity', 0));
            VR_headset_RR_matrix_visionpro(contrast_index, luminance_index) = binary_search_vec(bs_func, -find_sensitivity, [8 400], 20);
        end
    end
    writematrix(luminance_list_visionpro, 'VR_headset_results/luminance_list_visionpro.csv');
    writematrix(VR_headset_RR_matrix_visionpro, 'VR_headset_results/VR_headset_RR_matrix_visionpro.csv');
    
    width = 43;
    height = 29;
    luminance_list_hololens = logspace(log10(10),log10(10000),20);

    VR_headset_RR_matrix_hololens = zeros(length(contrast_list), length(luminance_list_hololens));
    for contrast_index = 1:length(contrast_list)
        contrast = contrast_list(contrast_index);
        find_sensitivity = 1 / contrast;
        for luminance_index = 1:length(luminance_list_hololens)
            luminance = luminance_list_hololens(luminance_index);
            bs_func = @(omega) - csf_elaTCSF_model.sensitivity_rect(struct('s_frequency', 0, 't_frequency', omega, 'orientation', 0, ...
                'luminance', luminance.*0.1, 'width', width, 'height', height, 'eccentricity', 0));
            VR_headset_RR_matrix_hololens(contrast_index, luminance_index) = binary_search_vec(bs_func, -find_sensitivity, [8 400], 20);
        end
    end
    writematrix(luminance_list_hololens, 'VR_headset_results/luminance_list_hololens.csv');
    writematrix(VR_headset_RR_matrix_hololens, 'VR_headset_results/VR_headset_RR_matrix_hololens.csv');

    width = 110;
    height = 96;
    luminance_list_metaquest = logspace(log10(10),log10(10000),20);

    VR_headset_RR_matrix_metaquest = zeros(length(contrast_list), length(luminance_list_metaquest));
    for contrast_index = 1:length(contrast_list)
        contrast = contrast_list(contrast_index);
        find_sensitivity = 1 / contrast;
        for luminance_index = 1:length(luminance_list_metaquest)
            luminance = luminance_list_metaquest(luminance_index);
            bs_func = @(omega) - csf_elaTCSF_model.sensitivity_rect(struct('s_frequency', 0, 't_frequency', omega, 'orientation', 0, ...
                'luminance', luminance.*0.1, 'width', width, 'height', height, 'eccentricity', 0));
            VR_headset_RR_matrix_metaquest(contrast_index, luminance_index) = binary_search_vec(bs_func, -find_sensitivity, [8 400], 20);
        end
    end
    writematrix(luminance_list_metaquest, 'VR_headset_results/luminance_list_metaquest.csv');
    writematrix(VR_headset_RR_matrix_metaquest, 'VR_headset_results/VR_headset_RR_matrix_metaquest.csv');
else
    luminance_list_visionpro = readmatrix('VR_headset_results/luminance_list_visionpro.csv');
    luminance_list_hololens = readmatrix('VR_headset_results/luminance_list_hololens.csv');
    luminance_list_metaquest = readmatrix('VR_headset_results/luminance_list_metaquest.csv');
    VR_headset_RR_matrix_visionpro = readmatrix('VR_headset_results/VR_headset_RR_matrix_visionpro.csv');
    VR_headset_RR_matrix_hololens = readmatrix('VR_headset_results/VR_headset_RR_matrix_hololens.csv');
    VR_headset_RR_matrix_metaquest = readmatrix('VR_headset_results/VR_headset_RR_matrix_metaquest.csv');
end

CFF_ticks = [40,50,60,70,80,90,100,110,120,130];
x_ticks_visio_pro = [10,100,1000,10000];
x_ticks_visio_hololens = [10,100,1000,10000]; %[10,100,500];
x_ticks_visio_metaquest = [10,100,1000,10000]; %[10,100];
if plot_need == 1
    ha = tight_subplot(1, 3, [.13 .05],[.12 .07],[.06 .02]);
    set(gcf, 'Position', [100, 100, 1000, 450]);
    axes(ha(1));
    % img = imread('E:\All_Conference_Papers\SIGGRAPH24/applevisionpro.png', 'png');
    % imshow(img);
    % hold on;
    x_fill = [luminance_list_visionpro, flip(luminance_list_visionpro)];
    y_fill = [VR_headset_RR_matrix_visionpro, ones(size(VR_headset_RR_matrix_visionpro)).*max(CFF_ticks)];
    fill(x_fill, y_fill, [1 0.8 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
    hold on;
    hh = [];
    plot(luminance_list_visionpro, VR_headset_RR_matrix_visionpro, 'Color', 'r', 'DisplayName', 'Apple Vision Pro', 'LineWidth', 2);
    hh(end+1) = plot(luminance_list_visionpro, ones(size(luminance_list_visionpro)).*90, 'Color', 'r', 'DisplayName', '90 Hz', 'LineStyle', '--');
    hh(end+1) = plot(luminance_list_visionpro, ones(size(luminance_list_visionpro)).*96, 'Color', 'r', 'DisplayName', '96 Hz', 'LineStyle', '--');
    hh(end+1) = plot(luminance_list_visionpro, ones(size(luminance_list_visionpro)).*100, 'Color', 'r', 'DisplayName', '100 Hz', 'LineStyle', '--');
    xline(peak_luminance_visio_pro, 'r--', 'LineWidth', 1);
    text(peak_luminance_visio_pro, min(CFF_ticks), ['peak lum: ' num2str(peak_luminance_visio_pro)], 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');
    x_center = 10^mean([log10(min(x_ticks_visio_pro)), log10(max(x_ticks_visio_pro))]);
    y_center = (min(CFF_ticks) + max(CFF_ticks))/2;
    text(x_center, 90, '90 Hz', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
    text(x_center, 96, '96 Hz', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
    text(x_center, 100, '100 Hz', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
    text(x_center, 110, 'FOV: 100$^{\circ}\times$100$^{\circ}$', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
        'Fontsize', 14, 'Interpreter', 'latex');
    % text(peak_luminance_visio_pro, y_center, ['Peak Lum. ' num2str(peak_luminance_visio_pro) 'cd/m^2'], 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
    xlabel('Luminance (cd/m^2)', FontSize=14);
    set(gca, 'XScale', 'log');
    ylabel('VR headset CFF (Hz)', FontSize=14);
    title('Apple Vision Pro (2024)', FontSize=14);
    % legend (hh, 'Location', 'best', 'Orientation', 'horizontal');
    xticks(x_ticks_visio_pro);
    xticklabels(x_ticks_visio_pro);
    xlim([min(x_ticks_visio_pro), max(x_ticks_visio_pro)]);
    yticks(CFF_ticks);
    yticklabels(CFF_ticks);
    ylim([min(CFF_ticks), max(CFF_ticks)]);
    

    axes(ha(2));
    x_fill = [luminance_list_hololens, flip(luminance_list_hololens)];
    y_fill = [VR_headset_RR_matrix_hololens, ones(size(VR_headset_RR_matrix_hololens)).*max(CFF_ticks)];
    fill(x_fill, y_fill, [0.8 1 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
    hold on;
    hh = [];
    plot(luminance_list_hololens, VR_headset_RR_matrix_hololens, 'Color', 'g', 'DisplayName', 'Microsoft HoloLens 2', 'LineWidth', 2);
    hh(end+1) = plot(luminance_list_hololens, ones(size(luminance_list_hololens)).*60, 'Color', 'g', 'DisplayName', '60 Hz', 'LineStyle', '--');
    xline(peak_luminance_holoLens_2, 'g--', 'LineWidth', 1);
    text(peak_luminance_holoLens_2, min(CFF_ticks), ['peak lum: ' num2str(peak_luminance_holoLens_2)], 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
    x_center = 10^mean([log10(min(x_ticks_visio_hololens)), log10(max(x_ticks_visio_hololens))]);
    y_center = (min(CFF_ticks) + max(CFF_ticks))/2;
    text(x_center, 60, '60 Hz', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
    text(x_center, 110, 'FOV: 43$^{\circ}\times$29$^{\circ}$', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
        'Fontsize', 14, 'Interpreter', 'latex');
    xlabel('Luminance (cd/m^2)', FontSize=14);
    set(gca, 'XScale', 'log');
    title('Microsoft HoloLens 2 (2019)', FontSize=14);
    % legend (hh, 'Location', 'best', 'Orientation', 'horizontal');
    xticks(x_ticks_visio_hololens);
    xticklabels(x_ticks_visio_hololens);
    xlim([min(x_ticks_visio_hololens), max(x_ticks_visio_hololens)]);
    yticks(CFF_ticks);
    yticklabels(CFF_ticks);
    ylim([min(CFF_ticks), max(CFF_ticks)]);
    


    axes(ha(3));
    x_fill = [luminance_list_metaquest, flip(luminance_list_metaquest)];
    y_fill = [VR_headset_RR_matrix_metaquest, ones(size(VR_headset_RR_matrix_metaquest)).*max(CFF_ticks)];
    fill(x_fill, y_fill, [0.8 0.8 1], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
    hold on;
    hh = [];
    plot(luminance_list_metaquest, VR_headset_RR_matrix_metaquest, 'Color', 'b', 'DisplayName', 'Meta Quest 3', 'LineWidth', 2);
    hh(end+1) = plot(luminance_list_metaquest, ones(size(luminance_list_metaquest)).*72, 'Color', 'b', 'DisplayName', '72 Hz', 'LineStyle', '--');
    hh(end+1) = plot(luminance_list_metaquest, ones(size(luminance_list_metaquest)).*80, 'Color', 'b', 'DisplayName', '80 Hz', 'LineStyle', '--');
    hh(end+1) = plot(luminance_list_metaquest, ones(size(luminance_list_metaquest)).*90, 'Color', 'b', 'DisplayName', '90 Hz', 'LineStyle', '--');
    hh(end+1) = plot(luminance_list_metaquest, ones(size(luminance_list_metaquest)).*120, 'Color', 'b', 'DisplayName', '120 Hz', 'LineStyle', '--');
    xline(peak_luminance_quest_3, 'b--', 'LineWidth', 1);
    text(peak_luminance_quest_3, min(CFF_ticks), ['peak lum: ' num2str(peak_luminance_quest_3)], 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');
    x_center = 10^mean([log10(min(x_ticks_visio_metaquest)), log10(max(x_ticks_visio_metaquest))]);
    y_center = (min(CFF_ticks) + max(CFF_ticks))/2;
    text(x_center, 72, '72 Hz', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
    text(x_center, 80, '80 Hz', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
    text(x_center, 90, '90 Hz', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
    text(x_center, 120, '120 Hz', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
    text(x_center, 110, 'FOV: 110$^{\circ}\times$96$^{\circ}$', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
        'Fontsize', 14, 'Interpreter', 'latex');
    xlabel('Luminance (cd/m^2)', FontSize=14);
    set(gca, 'XScale', 'log');
    title('Meta Quest 3 (2023)', FontSize=14);
    % legend (hh, 'Location', 'best', 'Orientation', 'horizontal');
    xticks(x_ticks_visio_metaquest);
    xticklabels(x_ticks_visio_metaquest);
    xlim([min(x_ticks_visio_metaquest), max(x_ticks_visio_metaquest)]);
    yticks(CFF_ticks);
    yticklabels(CFF_ticks);
    ylim([min(CFF_ticks), max(CFF_ticks)]);
end