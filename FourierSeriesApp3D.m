%% This is a matlab scriptthat loads up a user interfaceand it allows a 
% user to select different types of waveform: square wave, trapezoidal,
% triangular, or sawtooth as well as sine and cosine. The user is allowed 
% to specify the periodicity of the waveform. The tool will then allow the
% user to select the number of Fourier series coefficients to include in 
% the reconstruction of the waveform for its approximation (the number of 
% harmonics). The user can change the plot so that to see the fully 
% reconstructed waveform due to the superposition of all the component 
% sines and cosines that comprise it. The harmonics are plotted so you can 
% see the individual frequencies of each pure tone and get a better 
% understanding of the duality between time and frequency by looking at the
% Fourier series terms in a three dimensional plot allowing the user to 
% view components from orthogonal directions. The app provides user control
% of the 3d plot view azimuth and elevation to look at the component
% sine/cosines from different perspectives
% Written by:
%     Erick Oberstar (oberstar@engr.wisc.edu) 
%     Program Director - AI/ML, Electrification, and Mechatronics
%     University of Wisconsin - Madison
%     College of Engineering
%     Office of Interdiciplinary Professional Programs
%     Created on: 2026-02-26
%     Code generation help provided by perplexity.ai

classdef FourierSeriesApp3D < matlab.apps.AppBase
    % FourierSeriesApp3D
    % 3D view with X = harmonic index or frequency, Y = time, Z = amplitude.

    properties (Access = public)
        % Main figure window for the app
        UIFigure                 matlab.ui.Figure

        % Waveform type selection (Square, Triangle, Sawtooth, etc.)
        WaveformDropDown         matlab.ui.control.DropDown
        WaveformLabel            matlab.ui.control.Label

        % Slider to control number of Fourier harmonics in the approximation
        HarmonicsSlider          matlab.ui.control.Slider
        HarmonicsLabel           matlab.ui.control.Label

        % 3D axes showing original signal, individual harmonics, and sum
        TimeAxes3D               matlab.ui.control.UIAxes

        % Period T0 numeric input and label
        PeriodEditField          matlab.ui.control.NumericEditField
        PeriodLabel              matlab.ui.control.Label

        % Button to force recomputation when several controls are changed
        RecomputeButton          matlab.ui.control.Button

        % Checkbox to toggle plotting of individual harmonic components
        ShowHarmonicsCheckBox    matlab.ui.control.CheckBox

        % Dropdown to toggle X‑axis between harmonic index and frequency (Hz)
        YAxisModeDropDown        matlab.ui.control.DropDown   % actually X-axis mode now
        YAxisModeLabel           matlab.ui.control.Label

        % Display‑only fundamental frequency f0 corresponding to T0
        FrequencyLabel           matlab.ui.control.Label      % f0 label
        FrequencyEditField       matlab.ui.control.NumericEditField % f0 display

        % Sliders and labels to control 3D view azimuth and elevation
        AzimuthSlider            matlab.ui.control.Slider     % view azimuth
        AzimuthLabel             matlab.ui.control.Label
        ElevationSlider          matlab.ui.control.Slider     % view elevation
        ElevationLabel           matlab.ui.control.Label
    end

    properties (Access = private)
        % Internal time vector, sampling frequency, and fundamental period
        t   double
        fs  double
        T0  double
    end

    methods (Access = private)

        function setupSignalGrid(app)
            % Update internal period from GUI
            app.T0 = app.PeriodEditField.Value;
        
            % Fundamental frequency
            if app.T0 > 0
                f0 = 1/app.T0;
            else
                f0 = 0;
            end
        
            % Target samples per period (keeps waveforms well‑resolved
            % across a wide range of fundamental frequencies)
            samplesPerPeriod = 100;          % tweak 50–200 as you like
        
            if f0 > 0
                % Adaptive fs: keep at least 5000 Hz, but scale with f0
                app.fs = max(samplesPerPeriod * f0, 5000);
            else
                app.fs = 5000;
            end
        
            % Time vector over 3 periods of the waveform
            tmax   = 3*app.T0;
            app.t  = 0:1/app.fs:tmax;
        
            % Update f0 display
            if app.T0 > 0
                app.FrequencyEditField.Value = f0;
            else
                app.FrequencyEditField.Value = NaN;
            end
        end

        function x = generateWaveform(app)
            % Generate the chosen base waveform over the current time grid

            wtype = app.WaveformDropDown.Value;
            T0    = app.T0;
            t     = app.t;
            f0    = 1/T0;
            w0    = 2*pi*f0;

            switch wtype
                case 'Sine'
                    x = sin(w0*t);

                case 'Cosine'
                    x = cos(w0*t);

                case 'Square'
                    % 50% duty cycle square wave centered about zero[web:12]
                    x = square(w0*t);

                case 'Triangle'
                    % Triangle wave phase‑aligned with Fourier sine series
                    % implementation by shifting a symmetric sawtooth
                    x = sawtooth(w0*(t + T0/4), 0.5);

                case 'Sawtooth'
                    % Standard rising sawtooth on [0,T0]
                    x = sawtooth(w0*t);

                case 'Trapezoid'
                    % Trapezoid constructed by filtering a duty‑cycle square
                    duty = 0.4;             % flat‑top duty fraction
                    rise = 0.15;            % normalized rise time fraction
                    base = square(w0*t, duty*100);
                    winSamples = round(rise*T0*app.fs);
                    if winSamples < 1
                        winSamples = 1;
                    end
                    h  = ones(1, winSamples)/winSamples;
                    x  = conv(base, h, 'same');
                    x  = x / max(abs(x));   % normalize amplitude

                otherwise
                    % Fallback to a sine wave
                    x = sin(w0*t);
            end
        end

        function [xN, harmonics] = fourierApproxWithTerms(app, N)
            % Compute Fourier approximation xN using N harmonics and
            % return each individual harmonic term in 'harmonics'

            wtype = app.WaveformDropDown.Value;
            T0    = app.T0;
            t     = app.t;
            f0    = 1/T0;
            w0    = 2*pi*f0;
            harmonics = {};

            switch wtype
                case 'Square'
                    % Classic odd‑harmonic sine series for square wave[web:12]
                    xN = zeros(size(t));
                    for kk = 1:N
                        n    = 2*kk - 1;
                        term = (4/pi)*sin(n*w0*t)/n;
                        xN   = xN + term;
                        harmonics{kk} = term; %#ok<AGROW>
                    end

                case 'Triangle'
                    % Odd‑harmonic sine series with 1/n^2 decay
                    xN = zeros(size(t));
                    for kk = 1:N
                        n    = 2*kk - 1;
                        term = (8/(pi^2))*((-1)^(kk-1)*sin(n*w0*t)/n^2);
                        xN   = xN + term;
                        harmonics{kk} = term; %#ok<AGROW>
                    end

                case 'Sawtooth'
                    % Full‑harmonic sine series for sawtooth
                    xN = zeros(size(t));
                    for n = 1:N
                        term = -2/pi*sin(n*w0*t)/n;
                        xN   = xN + term;
                        harmonics{n} = term; %#ok<AGROW>
                    end

                case 'Sine'
                    % Single‑harmonic sine (no true series needed)
                    xN = sin(w0*t);
                    harmonics{1} = xN;

                case 'Cosine'
                    % Single‑harmonic cosine
                    xN = cos(w0*t);
                    harmonics{1} = xN;

                otherwise
                    % For trapezoid (and any other waveform falling here),
                    % derive harmonics from FFT and reconstruct each term[web:13]

                    % Sampled waveform over one grid
                    xFull = app.generateWaveform();
                    L     = numel(xFull);

                    % Complex Fourier coefficients (discrete‑time FS)
                    X     = fft(xFull)/L;

                    % Preallocate harmonic terms and total truncated spectrum
                    harmonics = cell(1,N);
                    XtruncAll = zeros(size(X));
                    XtruncAll(1) = X(1);          % DC component

                    for n = 1:N
                        % Build spectrum with only ±n harmonic present
                        Xn = zeros(size(X));

                        % Positive harmonic n at index n+1
                        if n+1 <= L
                            Xn(n+1) = X(n+1);
                        end

                        % Negative harmonic n at wrapped index (L-n+1)
                        idxNeg = mod(L-n+1-1, L) + 1;
                        Xn(idxNeg) = X(idxNeg);

                        % Time‑domain contribution of nth harmonic
                        harmonics{n} = real(ifft(Xn*L));

                        % Accumulate into truncated spectrum
                        XtruncAll = XtruncAll + Xn;
                    end

                    % Total truncated Fourier approximation
                    xN = real(ifft(XtruncAll*L));
            end
        end

        function applyView(app)
            % Apply azimuth/elevation settings from sliders to 3D axes
            az = app.AzimuthSlider.Value;
            el = app.ElevationSlider.Value;
            view(app.TimeAxes3D, az, el);
        end

        function plotTimeDomain3D(app, x, xN, harmonics)
            % Plot original signal, individual harmonics, and sum in 3D

            t  = app.t;
            T0 = app.T0;
            f0 = 1/app.T0;

            cla(app.TimeAxes3D);
            ax = app.TimeAxes3D;
            hold(ax, 'on');

            modeStr = app.YAxisModeDropDown.Value; % controls X-axis: index vs Hz

            % Original waveform at X = 0 (baseline sheet)
            x0 = zeros(size(t));
            plot3(ax, x0, t, x, 'b', 'LineWidth', 1.5);

            % Individual harmonics as vertical "sheets" over time
            if app.ShowHarmonicsCheckBox.Value && ~isempty(harmonics)
                dimColor = [0.6 0.6 0.6];
                for k = 1:numel(harmonics)
                    switch modeStr
                        case 'Harmonic index'
                            xVal = k;            % 1,2,3,...[web:71]
                        case 'Frequency (Hz)'
                            xVal = k*f0;         % physical frequency in Hz
                        otherwise
                            xVal = k;
                    end
                    xk = xVal*ones(size(t));
                    plot3(ax, xk, t, harmonics{k}, ...
                          'Color', dimColor, 'LineWidth', 0.8);
                end
            end

            % Summed Fourier approximation also at X = 0
            xSum = zeros(size(t));
            plot3(ax, xSum, t, xN, 'r--', 'LineWidth', 1.5);

            hold(ax, 'off');
            grid(ax, 'on');
            xlabel(ax, modeStr);
            ylabel(ax, 'Time (s)');
            zlabel(ax, 'Amplitude');
            title(ax, '3D Frequency / Time / Amplitude');

            % Use current slider values for 3D view orientation
            applyView(app);
        end

        function recomputeAndPlot(app)
            % Master update: rebuild time grid, regenerate waveform,
            % recompute Fourier series, and refresh plots
            setupSignalGrid(app);
            N  = round(app.HarmonicsSlider.Value);
            x  = app.generateWaveform();
            [xN, harmonics] = app.fourierApproxWithTerms(N);
            plotTimeDomain3D(app, x, xN, harmonics);
        end
    end

    % Callbacks
    methods (Access = private)

        function HarmonicsSliderValueChanged(app, src, event)
            %#ok<NASGU>  % src/event unused; callback only triggers recompute
            recomputeAndPlot(app);
        end

        function WaveformDropDownValueChanged(app, src, event)
            %#ok<NASGU>
            recomputeAndPlot(app);
        end

        function PeriodEditFieldValueChanged(app, src, event)
            %#ok<NASGU>
            recomputeAndPlot(app);
        end

        function RecomputeButtonPushed(app, src, event)
            %#ok<NASGU>
            recomputeAndPlot(app);
        end

        function ShowHarmonicsCheckBoxValueChanged(app, src, event)
            %#ok<NASGU>
            recomputeAndPlot(app);
        end

        function YAxisModeDropDownValueChanged(app, src, event)
            %#ok<NASGU>
            recomputeAndPlot(app);
        end

        function AzimuthSliderValueChanged(app, src, event)
            %#ok<NASGU>
            applyView(app);
        end

        function ElevationSliderValueChanged(app, src, event)
            %#ok<NASGU>
            applyView(app);
        end
    end

    % Component initialization
    methods (Access = private)

        function createComponents(app)
            % Create top‑level app window
            app.UIFigure = uifigure( ...
                'Name','Fourier Series Demo 3D', ...
                'WindowStyle','normal');
            app.UIFigure.Position = [100 100 900 600];

            % Waveform dropdown
            app.WaveformLabel = uilabel(app.UIFigure);
            app.WaveformLabel.Position = [25 560 80 22];
            app.WaveformLabel.Text = 'Waveform';
            app.WaveformDropDown = uidropdown(app.UIFigure);
            app.WaveformDropDown.Position = [110 560 140 22];
            app.WaveformDropDown.Items = ...
                {'Square','Sine','Cosine','Triangle','Sawtooth','Trapezoid'};
            app.WaveformDropDown.Value = 'Square';
            app.WaveformDropDown.ValueChangedFcn = ...
                @app.WaveformDropDownValueChanged;

            % Period T0 numeric input
            app.PeriodLabel = uilabel(app.UIFigure);
            app.PeriodLabel.Position = [280 560 80 22];
            app.PeriodLabel.Text = 'Period T0 (s)';
            app.PeriodEditField = uieditfield(app.UIFigure,'numeric');
            app.PeriodEditField.Position = [360 560 80 22];
            app.PeriodEditField.Value = 1;
            app.PeriodEditField.ValueChangedFcn = ...
                @app.PeriodEditFieldValueChanged;

            % Frequency display (f0) just below period
            app.FrequencyLabel = uilabel(app.UIFigure);
            app.FrequencyLabel.Position = [280 535 80 22];
            app.FrequencyLabel.Text = 'f0 (Hz)';
            app.FrequencyEditField = uieditfield(app.UIFigure,'numeric');
            app.FrequencyEditField.Position = [360 535 80 22];
            app.FrequencyEditField.Editable = 'off';

            % Harmonics slider (controls number of terms in series)
            app.HarmonicsLabel = uilabel(app.UIFigure);
            app.HarmonicsLabel.Position = [460 560 120 22];
            app.HarmonicsLabel.Text = '# Harmonics';
            app.HarmonicsSlider = uislider(app.UIFigure);
            app.HarmonicsSlider.Position = [560 570 300 3];
            app.HarmonicsSlider.Limits = [1 50];
            app.HarmonicsSlider.MajorTicks = [1 5 10 20 30 40 50];
            app.HarmonicsSlider.Value = 10;
            app.HarmonicsSlider.ValueChangedFcn = ...
                @app.HarmonicsSliderValueChanged;

            % Checkbox toggling display of individual harmonic traces
            app.ShowHarmonicsCheckBox = uicheckbox(app.UIFigure);
            app.ShowHarmonicsCheckBox.Position = [25 490 200 22];
            app.ShowHarmonicsCheckBox.Text = 'Show individual harmonics';
            app.ShowHarmonicsCheckBox.Value = 1;
            app.ShowHarmonicsCheckBox.ValueChangedFcn = ...
                @app.ShowHarmonicsCheckBoxValueChanged;[web:40]

            % X-axis mode dropdown (harmonic index vs physical frequency)
            app.YAxisModeLabel = uilabel(app.UIFigure);
            app.YAxisModeLabel.Position = [250 490 100 22];
            app.YAxisModeLabel.Text = 'X-axis mode';
            app.YAxisModeDropDown = uidropdown(app.UIFigure);
            app.YAxisModeDropDown.Position = [340 490 140 22];
            app.YAxisModeDropDown.Items = ...
                {'Harmonic index','Frequency (Hz)'};
            app.YAxisModeDropDown.Value = 'Frequency (Hz)';
            app.YAxisModeDropDown.ValueChangedFcn = ...
                @app.YAxisModeDropDownValueChanged;

            % Azimuth slider (3D view azimuth, -180..180 degrees)
            app.AzimuthLabel = uilabel(app.UIFigure);
            app.AzimuthLabel.Position = [25 460 80 22];
            app.AzimuthLabel.Text = 'Azimuth';
            app.AzimuthSlider = uislider(app.UIFigure);
            app.AzimuthSlider.Position = [110 470 300 3];
            app.AzimuthSlider.Limits = [-180 180];
            app.AzimuthSlider.Value = -45;
            app.AzimuthSlider.ValueChangedFcn = ...
                @app.AzimuthSliderValueChanged;

            % Elevation slider (3D view elevation, -90..90 degrees)
            app.ElevationLabel = uilabel(app.UIFigure);
            app.ElevationLabel.Position = [430 460 80 22];
            app.ElevationLabel.Text = 'Elevation';
            app.ElevationSlider = uislider(app.UIFigure);
            app.ElevationSlider.Position = [515 470 300 3];
            app.ElevationSlider.Limits = [-90 90];
            app.ElevationSlider.Value = 30;
            app.ElevationSlider.ValueChangedFcn = ...
                @app.ElevationSliderValueChanged;

            % Recompute button (explicit refresh)
            app.RecomputeButton = uibutton(app.UIFigure,'push');
            app.RecomputeButton.Position = [520 490 120 25];
            app.RecomputeButton.Text = 'Recompute';
            app.RecomputeButton.ButtonPushedFcn = ...
                @app.RecomputeButtonPushed;

            % 3D axes: X = frequency or harmonic index, Y = time, Z = amplitude
            app.TimeAxes3D = uiaxes(app.UIFigure);
            app.TimeAxes3D.Position = [50 50 800 360];
            title(app.TimeAxes3D, '3D Frequency / Time / Amplitude');
            xlabel(app.TimeAxes3D, 'Frequency (Hz)');
            ylabel(app.TimeAxes3D, 'Time (s)');
            zlabel(app.TimeAxes3D, 'Amplitude');
            grid(app.TimeAxes3D,'on');
        end
    end

    methods (Access = public)

        function app = FourierSeriesApp3D
            % App constructor: build UI, initialize grid, and draw plots
            set(0,'CurrentFigure',[]);
            createComponents(app);
            setupSignalGrid(app);
            recomputeAndPlot(app);
        end

        function delete(app)
            % Destructor: cleanly close the UIFigure when app is destroyed
            if isvalid(app.UIFigure)
                delete(app.UIFigure);
            end
        end
    end
end
