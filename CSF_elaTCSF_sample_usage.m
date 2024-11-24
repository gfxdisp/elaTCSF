csf_elaTCSF_model = CSF_elaTCSF();
csf_pars = struct('s_frequency', 0, 't_frequency', 10, 'orientation', 0, ...
                    'luminance', 3, 'area', 50.^2*pi, 'eccentricity', 10);
sensitivity_disk = csf_elaTCSF_model.sensitivity(csf_pars); % Sensitivity for Disk

csf_pars = struct('s_frequency', 0, 't_frequency', 10, 'orientation', 0, ...
                    'luminance', 3, 'width', 100, 'height', 100,'eccentricity', 10);
sensitivity_rect = csf_elaTCSF_model.sensitivity_rect(csf_pars); % Sensitivity for Rectangle
