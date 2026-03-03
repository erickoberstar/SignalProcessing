%% Fourier Series Demo — 2D (Simple) version
% GUI to demonstrate Fourier series approximations
% for several periodic waveforms, with time- and
% frequency-domain plots and a harmonics "knob".
% Written by:
%     Erick Oberstar (oberstar@engr.wisc.edu)
%     Program Director - AI/ML, Electrification, and Mechatronics
%     University of Wisconsin - Madison
%     College of Engineering
%     Office of Interdisciplinary Professional Programs
%     Created on: 2026-02-26
%     Code generation help provided by perplexity.ai
classdef FourierSeriesAppSimple < matlab.apps.AppBase

    properties (Access = public)
        UIFigure          matlab.ui.Figure
        WaveformDropDown  matlab.ui.control.DropDown
        WaveformLabel     matlab.ui.control.Label
        HarmonicsSlider   matlab.ui.control.Slider
        HarmonicsLabel    matlab.ui.control.Label
        TimeAxes          matlab.ui.control.UIAxes
        FreqAxes          matlab.ui.control.UIAxes
        PeriodEditField   matlab.ui.control.NumericEditField
        PeriodLabel       matlab.ui.control.Label
        RecomputeButton   matlab.ui.control.Button
    end

    properties (Access = private)
        t   double  % time vector
        fs  double  % sampling frequency
        T0  double  % fundamental period
    end

    methods (Access = private)

        function setupSignalGrid(app)
            app.T0 = app.PeriodEditField.Value;
            app.fs = 5000;
            tmax   = 3*app.T0;
            app.t  = 0:1/app.fs:tmax;
        end

        function x = generateWaveform(app)
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
                    x = square(w0*t);
                case 'Triangle'
                    x = sawtooth(w0*t, 0.5);
                case 'Sawtooth'
                    x = sawtooth(w0*t);
                case 'Trapezoid'
                    duty = 0.4;
                    rise = 0.15;
                    base = square(w0*t, duty*100);
                    winSamples = round(rise*T0*app.fs);
                    if winSamples < 1
                        winSamples = 1;
                    end
                    h = ones(1, winSamples)/winSamples;
                    x = conv(base, h, 'same');
                    x = x / max(abs(x));
                otherwise
                    x = sin(w0*t);
            end
        end

        function xN = fourierApprox(app, N)
            wtype = app.WaveformDropDown.Value;
            T0    = app.T0;
            t     = app.t;
            f0    = 1/T0;
            w0    = 2*pi*f0;
            switch wtype
                case 'Square'
                    xN = zeros(size(t));
                    for kk = 1:N
                        n  = 2*kk - 1;
                        xN = xN + sin(n*w0*t)/n;
                    end
                    xN = (4/pi)*xN;
                case 'Triangle'
                    xN = zeros(size(t));
                    for kk = 1:N
                        n  = 2*kk - 1;
                        xN = xN + ((-1)^((n-1)/2) * cos(n*w0*t) / n^2);
                    end
                    xN = (8/(pi^2)) * xN;
                case 'Sawtooth'
                    xN = zeros(size(t));
                    for n = 1:N
                        xN = xN + sin(n*w0*t)/n;
                    end
                    xN = -2/pi * xN;
                case 'Sine'
                    xN = sin(w0*t);
                case 'Cosine'
                    xN = cos(w0*t);
                otherwise
                    xFull = app.generateWaveform();
                    L     = numel(xFull);
                    X     = fft(xFull)/L;
                    Xtrunc      = zeros(size(X));
                    Xtrunc(1)   = X(1);
                    for n = 1:N
                        if n+1 <= L
                            Xtrunc(n+1) = X(n+1);
                        end
                        idxNeg = mod(L-n+1-1, L) + 1;
                        Xtrunc(idxNeg) = X(idxNeg);
                    end
                    xN = real(ifft(Xtrunc*L));
            end
        end

        function plotTimeDomain(app, x, xN)
            t  = app.t;
            T0 = app.T0;
            cla(app.TimeAxes);
            plot(app.TimeAxes, t, x, 'b', 'LineWidth', 1.5);
            hold(app.TimeAxes, 'on');
            plot(app.TimeAxes, t, xN, 'r--', 'LineWidth', 1.5);
            hold(app.TimeAxes, 'off');
            grid(app.TimeAxes, 'on');
            xlim(app.TimeAxes, [0, 3*T0]);
            xlabel(app.TimeAxes, 'Time (s)');
            ylabel(app.TimeAxes, 'Amplitude');
            title(app.TimeAxes, 'Time Domain: Original vs Fourier Approximation');
            legend(app.TimeAxes, {'Original', 'Approximation'}, 'Location', 'best');
        end

        function plotFreqDomain(app, N)
            % Plot analytic Fourier coefficients directly — no FFT leakage
            wtype = app.WaveformDropDown.Value;
            f0    = 1/app.T0;
            ax    = app.FreqAxes;
            nMax  = 2*N;  % plot out to 2N so zero even stems are visible

            cla(ax);
            hold(ax, 'on');

            switch wtype
                case 'Square'
                    for n = 1:nMax
                        freq = n*f0;
                        if mod(n,2) == 1          % odd: nonzero coefficient
                            cn = (4/pi)/n;
                            stem(ax, freq, cn, 'filled', 'MarkerSize', 5, ...
                                 'Color', [0.2 0.2 0.8], ...
                                 'MarkerFaceColor', [0.2 0.2 0.8]);
                        else                       % even: zero stem shown explicitly
                            stem(ax, freq, 0, 'filled', 'MarkerSize', 4, ...
                                 'Color', [0.75 0.75 0.85], ...
                                 'MarkerFaceColor', [0.75 0.75 0.85]);
                        end
                    end

                case 'Triangle'
                    for n = 1:nMax
                        freq = n*f0;
                        if mod(n,2) == 1
                            cn = (8/pi^2)/n^2;
                            stem(ax, freq, cn, 'filled', 'MarkerSize', 5, ...
                                 'Color', [0.2 0.2 0.8], ...
                                 'MarkerFaceColor', [0.2 0.2 0.8]);
                        else
                            stem(ax, freq, 0, 'filled', 'MarkerSize', 4, ...
                                 'Color', [0.75 0.75 0.85], ...
                                 'MarkerFaceColor', [0.75 0.75 0.85]);
                        end
                    end

                case 'Sawtooth'
                    for n = 1:N
                        cn = (2/pi)/n;
                        stem(ax, n*f0, cn, 'filled', 'MarkerSize', 5, ...
                             'Color', [0.2 0.2 0.8], ...
                             'MarkerFaceColor', [0.2 0.2 0.8]);
                    end

                case {'Sine','Cosine'}
                    stem(ax, f0, 1, 'filled', 'MarkerSize', 6, ...
                         'Color', [0.2 0.2 0.8], ...
                         'MarkerFaceColor', [0.2 0.2 0.8]);

                otherwise  % Trapezoid — use FFT (no closed-form series)
                    xN   = app.fourierApprox(N);
                    L    = numel(xN);
                    X    = fft(xN)/L;
                    k    = 0:L-1;
                    f    = k*app.fs/L;
                    fMax = (nMax+1)*f0;
                    idx  = f <= fMax;
                    stem(ax, f(idx), 2*abs(X(idx)), 'filled', 'MarkerSize', 4, ...
                         'Color', [0.2 0.2 0.8], ...
                         'MarkerFaceColor', [0.2 0.2 0.8]);
            end

            hold(ax, 'off');
            grid(ax, 'on');
            xlabel(ax, 'Frequency (Hz)');
            ylabel(ax, '|c_n|');
            title(ax, 'Fourier Coefficients (analytic)');
        end

        function recomputeAndPlot(app)
            setupSignalGrid(app);
            N  = round(app.HarmonicsSlider.Value);
            x  = app.generateWaveform();
            xN = app.fourierApprox(N);
            plotTimeDomain(app, x, xN);
            plotFreqDomain(app, N);   % <-- passes N, not xN
        end

    end

    % Callbacks
    methods (Access = private)

        function HarmonicsSliderValueChanged(app, src, event) %#ok
            recomputeAndPlot(app);
        end

        function WaveformDropDownValueChanged(app, src, event) %#ok
            recomputeAndPlot(app);
        end

        function PeriodEditFieldValueChanged(app, src, event) %#ok
            recomputeAndPlot(app);
        end

        function RecomputeButtonPushed(app, src, event) %#ok
            recomputeAndPlot(app);
        end

    end

    % Component initialization
    methods (Access = private)

        function createComponents(app)
            app.UIFigure = uifigure('Name', 'Fourier Series Demo');
            app.UIFigure.Position = [100 100 900 600];

            app.WaveformLabel = uilabel(app.UIFigure);
            app.WaveformLabel.Position = [25 560 80 22];
            app.WaveformLabel.Text = 'Waveform';

            app.WaveformDropDown = uidropdown(app.UIFigure);
            app.WaveformDropDown.Position = [110 560 140 22];
            app.WaveformDropDown.Items = ...
                {'Square','Sine','Cosine','Triangle','Sawtooth','Trapezoid'};
            app.WaveformDropDown.Value = 'Square';
            app.WaveformDropDown.ValueChangedFcn = @app.WaveformDropDownValueChanged;

            app.PeriodLabel = uilabel(app.UIFigure);
            app.PeriodLabel.Position = [280 560 80 22];
            app.PeriodLabel.Text = 'Period T0 (s)';

            app.PeriodEditField = uieditfield(app.UIFigure, 'numeric');
            app.PeriodEditField.Position = [360 560 80 22];
            app.PeriodEditField.Value = 1;
            app.PeriodEditField.ValueChangedFcn = @app.PeriodEditFieldValueChanged;

            app.HarmonicsLabel = uilabel(app.UIFigure);
            app.HarmonicsLabel.Position = [460 560 120 22];
            app.HarmonicsLabel.Text = '# Harmonics';

            app.HarmonicsSlider = uislider(app.UIFigure);
            app.HarmonicsSlider.Position = [560 570 300 3];
            app.HarmonicsSlider.Limits = [1 50];
            app.HarmonicsSlider.MajorTicks = [1 5 10 20 30 40 50];
            app.HarmonicsSlider.Value = 10;
            app.HarmonicsSlider.ValueChangedFcn = @app.HarmonicsSliderValueChanged;

            app.RecomputeButton = uibutton(app.UIFigure, 'push');
            app.RecomputeButton.Position = [25 520 120 25];
            app.RecomputeButton.Text = 'Recompute';
            app.RecomputeButton.ButtonPushedFcn = @app.RecomputeButtonPushed;

            app.TimeAxes = uiaxes(app.UIFigure);
            app.TimeAxes.Position = [50 280 800 220];
            title(app.TimeAxes, 'Time Domain');
            xlabel(app.TimeAxes, 'Time (s)');
            ylabel(app.TimeAxes, 'Amplitude');

            app.FreqAxes = uiaxes(app.UIFigure);
            app.FreqAxes.Position = [50 30 800 220];
            title(app.FreqAxes, 'Frequency Domain');
            xlabel(app.FreqAxes, 'Frequency (Hz)');
            ylabel(app.FreqAxes, '|c_n|');
        end

    end

    % App start and delete
    methods (Access = public)

        function app = FourierSeriesAppSimple
            createComponents(app);
            setupSignalGrid(app);
            recomputeAndPlot(app);
            registerApp(app, app.UIFigure);
        end

        function delete(app)
            delete(app.UIFigure);
        end

    end

end
