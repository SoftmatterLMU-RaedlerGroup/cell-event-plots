
%% Preparation for this script to run:
% load('datadir_test/corr_data_NP25_sts.mat')
% S = corrSTS(1);
% t_x(1) = eps;	% zero input not allowed for log-normal distribution
% distX = pdf(S.dist_x, t_x);
%
% Copyright © 2018 Daniel Woschée <daniel.woschee@physik.lmu.de>
% Faculty of Physics / Ludwig-Maximilians-Universität München
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License version 2
% as published by the Free Software Foundation.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program. If not, see <https://www.gnu.org/licenses/>.

%% Define fit models
logN = @(a, mu, sigma, x) a * ( exp(- (log(x)-mu).^2 /2 /sigma^2 ) ./x ) /sigma /sqrt(2*pi);
logN2 = @(a1, mu1, sigma1, a2, mu2, sigma2, x) logN(a1, mu1, sigma1, x) + logN(a2, mu2, sigma2, x);

logNormPdf = fittype(logN);
logNorm2Pdf = fittype(logN2);

logNormOpt = fitoptions(logNormPdf);
logNormOpt = fitoptions(logNormOpt, 'Lower', [0 eps eps], 'StartPoint', [1 1 1]);
logNorm2Opt = fitoptions(logNorm2Pdf);
logNorm2Opt = fitoptions(logNorm2Opt, 'Lower', [0 eps eps 0 eps eps], 'StartPoint', [.5 0 1 .5 1 1]);

%% Perform fit
foX = fit(t_x(:), distX(:), logNormPdf, logNormOpt);
foX2 = fit(t_x(:), distX(:), logNorm2Pdf, logNorm2Opt);

foX_raw = fit(t_x(:), S.counts_x(:), logNormPdf, logNormOpt);
foX2_raw = fit(t_x(:), S.counts_x(:), logNorm2Pdf, logNorm2Opt);

%% Get fit parameters
coeffs = num2cell(coeffvalues(foX2_raw));
[a1, mu1, sigma1, a2, mu2, sigma2] = coeffs{:};

%% Plot
fh = figure;
ax = axes(fh);
plot(ax, t_x, S.counts_x)
hold on
plot(ax, t_x, foX_raw(t_x))
plot(ax, t_x, foX2_raw(t_x))

plot(ax, t_x, logN(a1, mu1, sigma1, t_x))
plot(ax, t_x, logN(a2, mu2, sigma2, t_x))

title(ax, 'Event time correlations for tmrm-ros (STS)')
xlabel(ax, 't(event, tmrm) [h]')
ylabel(ax, 'Counts [#]')

%% Output
fh.PaperUnits = 'normalized';
fh.PaperPositionMode = 'manual';
fh.PaperPosition = [0 0 1 1];
fh.PaperType = 'A4';
fh.PaperOrientation = 'landscape';
print(fh, fullfile('out', 'lognormal_STS_tmrm-ros.pdf'), '-dpdf')