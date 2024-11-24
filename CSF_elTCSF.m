classdef CSF_elTCSF < CSF_base
    properties( Constant )
    end

    properties
        use_gpu = true;
    end

    methods

        function obj = CSF_elTCSF(  )
            obj.par = obj.get_default_par();
        end

        function name = short_name( obj )
            name = 'elTCSF';
        end

        function name = full_name( obj )
            name = 'elTCSF';
        end

        function S = TCSF(obj, omega, eccentricity, luminance)
            TCSF_n1 = 15;
            TCSF_n2 = 16;
            TCSF_xi = obj.par.TCSF_xi;
            TCSF_tau = obj.par.TCSF_tau;
            TCSF_kappa = obj.par.TCSF_kappa;
            TCSF_zeta = obj.par.TCSF_zeta;

            tcsf_lum_k1 = obj.par.tcsf_lum_k1;
            tcsf_ecc_k1 = obj.par.tcsf_ecc_k1;
            ecc_peak_f = obj.par.ecc_peak_f;
            omega = (omega - ecc_peak_f) ./ ((1+tcsf_ecc_k1.*eccentricity)) + ecc_peak_f;
            lum_peak_f = 0;
            omega = (omega - lum_peak_f) ./ (obj.par.tcsf_lum_b1+tcsf_lum_k1 * log10(luminance)) + lum_peak_f;
            S = abs(TCSF_xi * ((1 + 2 * 1i * pi * omega * TCSF_tau).^(-TCSF_n1) - TCSF_zeta * (1 + 2 * 1i * pi * omega * TCSF_kappa * TCSF_tau).^(-TCSF_n2)));
        end

        function S_ecc_factor = S_ecc(obj, eccentricity)
            ecc_k1 = obj.par.ecc_k1;
            S_ecc_factor = 10 .^ (-ecc_k1.*eccentricity);
        end

        function S_luminance_factor = S_lum(obj, luminance, eccentricity)
            S_luminance_factor = obj.par.lum_k(1) .* ((1+obj.par.lum_k(2)./luminance) .^ (-obj.par.lum_k(3)));
        end

        function S = sensitivity( obj, csf_pars )
            csf_pars = obj.test_complete_params(csf_pars, { 'luminance', 'ge_sigma' } );
            omega = csf_pars.t_frequency;
            s_frequency = csf_pars.s_frequency;
            eccentricity = csf_pars.eccentricity;
            luminance = csf_pars.luminance;
            S = s_frequency * 0 + obj.S_ecc(eccentricity) .* obj.S_lum(luminance, eccentricity) .* obj.TCSF(omega, eccentricity, luminance);
        end

        function S = sensitivity_edge( obj, csf_pars )
            csf_pars = obj.test_complete_params(csf_pars, { 'luminance', 'ge_sigma' } );
            omega = csf_pars.t_frequency;
            luminance = csf_pars.luminance;
            eccentricity = csf_pars.eccentricity;
            S = obj.S_ecc(eccentricity) .* obj.S_lum(luminance, eccentricity) .* obj.TCSF(omega, eccentricity, luminance);
            S = permute(S, circshift(1:numel(size(S)), -1));
        end

        function S = sensitivity_plot( obj, omega, luminance, eccentricity)
            S = obj.S_ecc(eccentricity) .* obj.S_lum(luminance, eccentricity) .* obj.TCSF(omega, eccentricity, luminance);
        end

        function pd = get_plot_description( obj )
            pd = struct();
            pp = 1;
            pd(pp).title = 'elTCSF11-Sensitivity-Temporal Frequency - different eccentricity';
            pd(pp).id = 'S_tf_ecc';
            pp = pp+1;

            pd(pp).title = 'elTCSF11-Sensitivity-Temporal Frequency - different luminance';
            pd(pp).id = 'S_tf_lum';
            pp = pp+1;

            pd(pp).title = 'elTCSF11-Sensitivity-Luminance';
            pd(pp).id = 'S_lum';
            pp = pp+1;

            pd(pp).title = 'elTCSF11-Sensitivity-Eccentricity';
            pd(pp).id = 'S_ecc';
            pp = pp+1;

            pd(pp).title = 'Ferry Porter Law';
            pd(pp).id = 'ferry_porter';
            pp = pp+1;
        end

        function plot_mechanism( obj, plt_id )
            switch( plt_id )
                case 'S_tf_ecc'
                    clf;
                    html_change_figure_print_size( gcf, 10, 10 );
                    omega = linspace( 0, 100 , 100)';
                    luminance = 3;
                    eccentricity_list = [0,10,20,30,40,50,60];
                    hh = cell(length(eccentricity_list), 1);
                    for eccentricity_index = 1:length(eccentricity_list)
                        eccentricity = eccentricity_list(eccentricity_index);
                        S_response = obj.sensitivity_plot(omega, luminance, eccentricity);
                        hh{eccentricity_index} = plot( omega, S_response, 'DisplayName', ['elTCSF - lum = 3 - ecc = ', num2str(eccentricity)]);
                        hold on;
                    end
                    xlabel( 'Temp. freq. [Hz]' );
                    ylabel( 'Sensitivity' );
                    set(gca, 'YScale', 'log');
                    legend([hh{:}], 'Location', 'Best');
                    grid on;

                case 'S_tf_lum'
                    clf;
                    html_change_figure_print_size( gcf, 10, 10 );
                    omega = linspace( 0, 100 , 100)';
                    luminance_list = [0.1,1,10,100];
                    eccentricity = 0;
                    hh = cell(length(luminance_list), 1);
                    for luminance_index = 1:length(luminance_list)
                        luminance = luminance_list(luminance_index);
                        S_response = obj.sensitivity_plot(omega, luminance, eccentricity);
                        hh{luminance_index} = plot( omega, S_response, 'DisplayName', ['elTCSF - ecc = 0 - lum = ', num2str(luminance)]);
                        hold on;
                    end
                    xlabel( 'Temp. freq. [Hz]' );
                    ylabel( 'Sensitivity' );
                    set(gca, 'YScale', 'log');
                    legend([hh{:}], 'Location', 'Best');
                    grid on;

                case 'S_lum'
                    clf;
                    html_change_figure_print_size( gcf, 10, 10 );
                    omega = 10;
                    luminance = logspace(log10(0.01), log10(100), 100)';
                    eccentricity_list = [0, 10, 20, 30, 40];
                    hh = [];
                    for eccentricity_index = 1:length(eccentricity_list)
                        eccentricity = eccentricity_list(eccentricity_index);
                        S_response = obj.sensitivity_plot(omega, luminance, eccentricity);
                        hh(end+1) = plot( luminance, S_response, 'DisplayName', ['elTCSF - omega = 10 - ecc = ' num2str(eccentricity)]);
                        hold on;
                    end

                    xlabel( 'Luminance (cd/m^2)' );
                    ylabel( 'Sensitivity' );
                    set(gca, 'YScale', 'log');
                    set(gca, 'XScale', 'log');
                    legend( hh, 'Location', 'Best' );
                    grid on;

                case 'S_ecc'
                    clf;
                    html_change_figure_print_size(gcf, 10, 10);
                    omega = 10;
                    luminance_list = [0.03, 0.3, 3, 30];
                    hh = cell(length(luminance_list), 1); % Corrected hh initialization
                    eccentricity = linspace(0, 90, 90)';
                    for luminance_index = 1:length(luminance_list)
                        luminance = luminance_list(luminance_index);
                        S_response = obj.sensitivity_plot(omega, luminance, eccentricity); % Corrected luminance variable name
                        hh{luminance_index} = plot(eccentricity, S_response, 'DisplayName', ['elTCSF - omega = 10 - lum = ', num2str(luminance)]);
                        hold on; % Added to keep multiple plots
                    end
                    hold off; % Added to release hold after plotting
                    xlabel('Eccentricity (degree)');
                    ylabel('Sensitivity');
                    set(gca, 'YScale', 'log');
                    legend([hh{:}], 'Location', 'Best'); % Corrected legend call
                    grid on;

                case 'ferry_porter'
                    clf;
                    html_change_figure_print_size( gcf, 10, 10 );
                    luminance = logspace( 0, 4, 200 )';
                    ECCs = [0 10 20 30];
                    hh = [];
                    for kk=1:length(ECCs)
                        bs_func = @(omega) -obj.find_cff(omega, luminance, ECCs(kk));
                        CFF = binary_search_vec(bs_func, -ones(numel(luminance),1), [10 160], 20);
                        hh(kk) = plot( luminance, CFF, 'DisplayName', sprintf( 'ecc=%g deg', ECCs(kk)) );
                        hold on
                    end
                    legend( hh, 'Location', 'best' );
                    set_axis_tick_label('x', 'luminance', luminance );
                    ylabel( 'CFF [Hz]' );

                otherwise
                    error( 'Wrong plt_id' );

            end
        end

        function S = find_cff( obj, omega, L, eccentricity)
            csf_pars = struct( 's_frequency', 0, 't_frequency', omega, 'orientation', 0, 'luminance', L, 'area', 1, 'eccentricity', eccentricity );
            S = obj.sensitivity( csf_pars );
        end

        function obj = set_pars( obj, pars_vector )
            obj = obj.set_pars@CSF_base(pars_vector);
            obj = obj.update_parameters();
        end

        function obj = update_parameters (obj)
            obj.par = CSF_base.update_struct( obj.par, obj.par );
        end

        function print( obj, fh )
            % Print the model parameters in a format ready to be pasted into
            % get_default_par()
            fprintf(fh, 'Parameters for TCSF component:\n');
            obj.print@CSF_base(fh)

            % Printed formatted parameters for the component classes
            % fprintf(fh, 'Parameters for TCSF component:\n');
            % obj.print(fh);
            % fprintf(fh, '\n');
        end
    end

    methods ( Static )

        function p = get_default_par()
            p = CSF_base.get_dataset_par();
            %% Parameters in the SIGGRAPH ASIA 2024
            p.lum_k = [ 1.76801 1.62402 0.533781 ];
            p.ecc_k1 = 0.0330933;
            p.tcsf_ecc_k1 = 0.0341811;
            p.tcsf_lum_k1 = 0.222269;
            p.tcsf_lum_b1 = 0.6678;
            p.TCSF_xi = 154.133;
            p.TCSF_tau = 0.00292069;
            p.TCSF_kappa = 2.12547;
            p.TCSF_zeta = 0.721095;
            % p.TCSF_n1 = 15;
            % p.TCSF_n2 = 16;
            p.ecc_peak_f = -2;
        end
    end
end
