classdef CSF_elaTCSF < CSF_base
    properties( Constant )
    end

    properties
        use_gpu = true;
        elTCSF = [];
    end

    methods

        function obj = CSF_elaTCSF(  )
            obj.par = obj.get_default_par();
            obj.elTCSF = CSF_elTCSF();
        end

        function name = short_name( obj )
            name = 'elaTCSF';
        end

        function name = full_name( obj )
            name = 'elaTCSF';
        end

        function S = sensitivity(obj, csf_pars )
            csf_pars = obj.test_complete_params(csf_pars, { 'luminance', 'ge_sigma' } );
            s_frequency = csf_pars.s_frequency;
            t_frequency = csf_pars.t_frequency;
            luminance = csf_pars.luminance;
            radius = csf_pars.ge_sigma;
            eccentricity = csf_pars.eccentricity;

            vector_length = max([length(s_frequency), length(t_frequency), length(luminance), length(radius), length(eccentricity)]);
            csf_pars_after_check = obj.check_and_fill_csf_pars({s_frequency, t_frequency, luminance, radius, eccentricity}, vector_length);
            s_frequency = csf_pars_after_check{1};
            t_frequency = csf_pars_after_check{2};
            luminance = csf_pars_after_check{3};
            radius = csf_pars_after_check{4};
            eccentricity = csf_pars_after_check{5};
            S = zeros(vector_length,1);
            for index = 1:vector_length
                S(index) = obj.Energy_S_ecc(obj.elTCSF, s_frequency(index), t_frequency(index), luminance(index), ...
                    radius(index), obj.par.E_thr, obj.par.beta, eccentricity(index));
            end
        end

        function S = sensitivity_rect(obj, csf_pars)
            csf_pars = obj.test_complete_params(csf_pars, { 'luminance'} );
            s_frequency = csf_pars.s_frequency;
            t_frequency = csf_pars.t_frequency;
            luminance = csf_pars.luminance;
            width = csf_pars.width;
            height = csf_pars.height;
            eccentricity = csf_pars.eccentricity;

            vector_length = max([length(s_frequency), length(t_frequency), length(luminance), length(width), length(height), length(eccentricity)]);
            csf_pars_after_check = obj.check_and_fill_csf_pars({s_frequency, t_frequency, luminance, width, height, ...
                eccentricity}, vector_length);
            s_frequency = csf_pars_after_check{1};
            t_frequency = csf_pars_after_check{2};
            luminance = csf_pars_after_check{3};
            width = csf_pars_after_check{4};
            height = csf_pars_after_check{5};
            eccentricity = csf_pars_after_check{6};
            S = zeros(vector_length,1);
            for index = 1:vector_length
                S(index) = obj.Energy_S_ecc_rect(obj.elTCSF, s_frequency(index), t_frequency(index), luminance(index), ...
                    width(index), height(index), obj.par.E_thr, obj.par.beta, eccentricity(index));
            end
        end

        function S = sensitivity_edge(obj, csf_pars )
            csf_pars = obj.test_complete_params(csf_pars, { 'luminance', 'lms_bkg', 'lms_delta', 'ge_sigma' } );
            t_frequency = csf_pars.t_frequency';
            s_frequency = zeros(size(t_frequency));
            luminance = csf_pars.luminance';
            radius = csf_pars.ge_sigma';
            eccentricity = csf_pars.eccentricity';

            vector_length = max([length(s_frequency), length(t_frequency), length(luminance), length(radius), length(eccentricity)]);
            csf_pars_after_check = obj.check_and_fill_csf_pars({s_frequency, t_frequency, luminance, radius, eccentricity}, vector_length);
            s_frequency = csf_pars_after_check{1};
            t_frequency = csf_pars_after_check{2};
            luminance = csf_pars_after_check{3};
            radius = csf_pars_after_check{4};
            eccentricity = csf_pars_after_check{5};
            S = zeros(vector_length,1);
            for index = 1:vector_length
                S(index) = obj.Energy_S_ecc(obj.elTCSF, s_frequency(index), t_frequency(index), luminance(index), ...
                    radius(index), obj.par.E_thr, obj.par.beta, eccentricity(index));
            end
        end

        function value = S_CSF(obj, csf_model, s_frequency, t_frequency, luminance, area, eccentricity)
            csf_pars = struct('s_frequency', s_frequency, 't_frequency', t_frequency, 'orientation', 0, ...
                'luminance', luminance, 'area', area, 'eccentricity', eccentricity);
            value = csf_model.sensitivity(csf_pars);
        end

        function energy_s = Energy_S_ecc_rect(obj, csf_model, s_frequency, t_frequency, luminance, width, height, E_thr, beta, eccentricity)
            S_ecc = @(degree_x, degree_y) obj.S_CSF(csf_model, s_frequency, t_frequency, luminance, 1, (degree_x.^2 + degree_y.^2).^0.5).^beta;
            intergration_value = integral2(S_ecc, eccentricity - width ./ 2, eccentricity + width ./ 2, -height ./ 2, height ./ 2);
            contrast = (E_thr ./ intergration_value).^(1/beta);
            energy_s = 1 ./ contrast;
        end

        function energy_s = Energy_S_ecc(obj, csf_model, s_frequency, t_frequency, luminance, radius, E_thr, beta, eccentricity)
            S_ecc = @(r,theta) obj.S_CSF(csf_model, s_frequency, t_frequency, luminance, 1, (r.^2 + eccentricity.^2 + 2.*eccentricity.*r.*cos(theta)).^0.5).^beta.*r;
            intergration_value = integral2(S_ecc, 0, radius, 0, 2*pi);
            contrast = (E_thr ./ intergration_value).^(1/beta);
            energy_s = 1 ./ contrast;
        end

        function csf_pars = check_and_fill_csf_pars(obj, csf_pars, v_length)
            for csf_par_i = 1:numel(csf_pars)
                csf_par = csf_pars{csf_par_i};
                if length(csf_par) == v_length
                    continue;
                elseif length(csf_par) == 1 && v_length ~= 1
                    csf_pars{csf_par_i} = repmat(csf_par, v_length, 1);
                else
                    disp('WRONG!');
                end
            end
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

            % pd(pp).title = 'Ferry Porter Law';
            % pd(pp).id = 'ferry_porter';
            % pp = pp+1;
        end

        function plot_mechanism( obj, plt_id )
            switch( plt_id )
                case 'S_tf_ecc'
                    clf;
                    html_change_figure_print_size( gcf, 10, 10 );
                    omega = linspace( 0, 100 , 100)';
                    luminance = 3;
                    area = 1;
                    eccentricity_list = [0,10,20,30,40,50,60];
                    hh = cell(length(eccentricity_list), 1);
                    for eccentricity_index = 1:length(eccentricity_list)
                        eccentricity = eccentricity_list(eccentricity_index);
                        csf_pars = struct('s_frequency', 0, 't_frequency', omega, 'orientation', 0, ...
                            'luminance', luminance, 'area', area, 'eccentricity', eccentricity);
                        S_response = obj.sensitivity(csf_pars);
                        hh{eccentricity_index} = plot( omega, S_response, 'DisplayName', ['elaTCSF - lum = 3; area = 1; ecc = ', num2str(eccentricity)]);
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
                    area = 1;
                    eccentricity = 0;
                    hh = cell(length(luminance_list), 1);
                    for luminance_index = 1:length(luminance_list)
                        luminance = luminance_list(luminance_index);
                        csf_pars = struct('s_frequency', 0, 't_frequency', omega, 'orientation', 0, ...
                            'luminance', luminance, 'area', area, 'eccentricity', eccentricity);
                        S_response = obj.sensitivity(csf_pars);
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
                    area = 1;
                    eccentricity = 0;
                    csf_pars = struct('s_frequency', 0, 't_frequency', omega, 'orientation', 0, ...
                        'luminance', luminance, 'area', area, 'eccentricity', eccentricity);
                    S_response = obj.sensitivity(csf_pars);
                    hh = plot( luminance, S_response, 'DisplayName', 'elTCSF - omega = 10 - ecc = 0');
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
                    area = 1;
                    eccentricity = linspace(0, 90, 90)';
                    for luminance_index = 1:length(luminance_list)
                        luminance = luminance_list(luminance_index);
                        csf_pars = struct('s_frequency', 0, 't_frequency', omega, 'orientation', 0, ...
                            'luminance', luminance, 'area', area, 'eccentricity', eccentricity);
                        S_response = obj.sensitivity(csf_pars);
                        hh{luminance_index} = plot(eccentricity, S_response, 'DisplayName', ['elTCSF - omega = 10 - lum = ', num2str(luminance)]);
                        hold on; % Added to keep multiple plots
                    end
                    hold off; % Added to release hold after plotting
                    xlabel('Eccentricity (degree)');
                    ylabel('Sensitivity');
                    set(gca, 'YScale', 'log');
                    legend([hh{:}], 'Location', 'Best'); % Corrected legend call
                    grid on;

                    % case 'ferry_porter'
                    %     clf;
                    %     html_change_figure_print_size( gcf, 10, 10 );
                    %     luminance = logspace( 0, 4, 200 )';
                    %     ECCs = [0 10 20 30];
                    %     hh = [];
                    %     for kk=1:length(ECCs)
                    %         bs_func = @(omega) -obj.find_cff(omega, luminance, ECCs(kk));
                    %         CFF = binary_search_vec(bs_func, -ones(numel(luminance),1), [10 160], 20);
                    %         hh(kk) = plot( luminance, CFF, 'DisplayName', sprintf( 'ecc=%g deg', ECCs(kk)) );
                    %         hold on
                    %     end
                    %     legend( hh, 'Location', 'best' );
                    %     set_axis_tick_label('x', 'luminance', luminance );
                    %     ylabel( 'CFF [Hz]' );

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

        % Copy parameters to from this object to individual CSF components
        function obj = update_parameters (obj)
            obj.par = CSF_base.update_struct( obj.par, obj.par );
            obj.elTCSF.par =  CSF_base.update_struct( obj.par.elTCSF,  obj.elTCSF.par);
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
            p.E_thr = 6.52801;
            p.beta = 3.80022;
            p.elTCSF.lum_k = [ 1.76801 1.62402 0.533781 ];
            p.elTCSF.ecc_k1 = 0.0330933;
            p.elTCSF.tcsf_ecc_k1 = 0.0341811;
            p.elTCSF.tcsf_lum_k1 = 0.222269;
            p.elTCSF.tcsf_lum_b1 = 0.6678;
            p.elTCSF.TCSF_xi = 154.133;
            p.elTCSF.TCSF_tau = 0.00292069;
            p.elTCSF.TCSF_kappa = 2.12547;
            p.elTCSF.TCSF_zeta = 0.721095;
        end
    end
end
