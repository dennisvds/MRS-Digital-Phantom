function [out, out_components, names] = create_lipids(n, sw, Bfield, linewidth, TE, seq, cf)
%CREATE_LIPIDS Create lipid spectra for different sequences
%
% Inputs:
%   n         - Number of points
%   sw        - Spectral width (Hz)
%   Bfield    - Magnetic field strength (Tesla)
%   linewidth - Linewidth (Hz)
%   TE        - Echo time (ms)
%   seq       - Sequence type ('PRESS', 'sLASER', 'LASER', 'STEAM', etc.)
%   cf        - Center frequency (Hz)
%
% Outputs:
%   out            - Summed lipid spectrum
%   out_components - Cell array of individual lipid spectra
%   names          - Cell array of lipid names (extracted from sys)

    % Define FA chain lengths and DB positions for spin systems
    % Simulate oleic acid (18:1 n-9)
    sys_oleic = create_sys(18, [9]);
    
    % Simulate palmitic acid (16:0)
    sys_palmitic = create_sys(16, []);
    
    % Simulate linoleic acid (18:2 n-6)
    sys_linoleic = create_sys(18, [6,9]);
    
    % Simulate palmitoleic acid (16:1 n-7)
    sys_palmitoleic = create_sys(16, [7]);

    % Store all systems in a list
    syss = {sys_oleic, sys_palmitic, sys_linoleic, sys_palmitoleic};

    % Lipid component names
    names = {'oleic','palmitic','linoleic','palmitoleic'};

    % Calculate center frequency in ppm
    gamma = 42.57747892; % MHz/T for 1H
    cf_ppm = cf / (Bfield * gamma); % cf in ppm

    % Initialize outputs
    out = [];
    out_components = cell(1, length(syss));

    % Set sequence to 'PRESS' as default (other sequenced not used in current implementation)
    seq = 'PRESS'; 

    % Loop over each component
    for k = 1:length(syss)
        sys_k = syss{k};  % Take single component

        % Simulate depending on the sequence
        switch lower(seq)
            case 'press'
                tau1 = 11.7;
                tau2 = TE-tau1;
                out_k = sim_press(n, sw, Bfield, linewidth, sys_k, tau1, tau2);

            case 'slaser'
                % Load adiabatic RF pulse
                rfPulse = io_loadRFwaveform('sampleAFPpulse_HS2_R15.RF', 'inv');
                refTp = 3.5; % RF pulse duration [ms]
                flipAngle = 180; % Flip angle [deg]

                % Spatial grid parameters
                thkX = 2; thkY = 2; % cm
                fovX = 3; fovY = 3; % cm
                nX = 16; nY = 16;
                x = linspace(-fovX/2, fovX/2, nX);
                y = linspace(-fovY/2, fovY/2, nY);

                % Resample pulse for speed
                rfPulse = rf_resample(rfPulse, 100);

                % Calculate gradients
                if ~rfPulse.isGM
                    Gx = (rfPulse.tbw/(refTp/1000))/(gamma*thkX/10000); %[G/cm]
                    Gy = (rfPulse.tbw/(refTp/1000))/(gamma*thkY/10000); %[G/cm]
                else
                    Gx = (rfPulse.tthk/(refTp/1000))/thkX;
                    Gy = (rfPulse.tthk/(refTp/1000))/thkY;
                end

                % Simulate semi-LASER for single sys_k
                out_k = [];
                for X = 1:length(x)
                    for Y = 1:length(y)
                        temp = sim_semiLASER_shaped(n, sw, Bfield, linewidth, sys_k, TE, ...
                            rfPulse, refTp, x(X), y(Y), Gx, Gy, flipAngle, cf_ppm);
                        out_k = op_addScans(out_k, temp);
                    end
                end

                % Scaling
                numSims = nX * nY;
                out_k = op_ampScale(out_k, 1/numSims);
                voxRatio = (thkX*thkY)/(fovX*fovY);
                out_k = op_ampScale(out_k, 1/voxRatio);

            case 'laser'
                out_k = sim_laser(n, sw, Bfield, linewidth, sys_k, TE);

            otherwise
                error('Unsupported sequence type: %s', seq);
        end

        % Store individual component
        out_components{k} = out_k;

        % Add to total output
        if isempty(out)
            out = out_k;
        else
            out = op_addScans(out, out_k);
        end
    end

    % % Plot the results
    % figure;
    % hold on;
    % for k = 1:length(out_components)
    %     % Extract specs from the struct
    %     specs = out_components{k}.specs;  % Assuming 'specs' is a field in the struct
    %     ppm = out_components{k}.ppm;
    %     plot(ppm, real(specs), 'DisplayName', names{k}, 'LineWidth', 1.2);  % Assuming 'f' is frequency and 'data' is the spectrum
    % end
    % hold off;
    % 
    % % Add labels and legend
    % xlabel('Chemical Shift (ppm)');
    % ylabel('Amplitude');
    % legend('show');
    % title('Lipid Spectra');
    % 
    % figure;
    % spec = out.specs;
    % ppm = flip(out.ppm);
    % plot(ppm, real(spec), 'LineWidth', 1.5)
    % % Add labels and legend
    % xlabel('Chemical Shift (ppm)');
    % ylabel('Amplitude');
    % title('Summed Lipid Spectra');
end
