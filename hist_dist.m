function varargout = hist_dist(sample, v_ev, is_normalized)
%HIST_DIST returns properties of a lognormal distribution of a sample
%
% Syntax
%   pdf_max = HIST_DIST(sample)
%   [pdf_max, err_max] = HIST_DIST(sample)
%   [pdf_max, err_max, peak] = HIST_DIST(sample)
%   [pdf_max, err_max, mu, sigma] = HIST_DIST(sample)
%   [pdf_max, err_max, peak, mu, sigma] = HIST_DIST(sample)
%   [pdf_max, err_max, mu, err_mu, sigma, err_sigma] = HIST_DIST(sample)
%   [pdf_max, err_max, peak, mu, err_mu, sigma, err_sigma] = HIST_DIST(sample)
%
%   [pdf_max, pdf] = HIST_DIST(sample, v_ev)
%   [pdf_max, pdf] = HIST_DIST(sample, v_ev, is_normalized)
%
%   [pdf_max, err_max, pdf] = HIST_DIST(sample, v_ev, ...)
%   [pdf_max, err_max, pdf, peak] = HIST_DIST(sample, v_ev, ...)
%   [pdf_max, err_max, pdf, mu, sigma] = HIST_DIST(sample, v_ev, ...)
%   [pdf_max, err_max, pdf, peak, mu, sigma] = HIST_DIST(sample, v_ev, ...)
%   [pdf_max, err_max, pdf, mu, err_mu, sigma, err_sigma] = HIST_DIST(sample, v_ev, ...)
%   [pdf_max, err_max, pdf, peak, mu, err_mu, sigma, err_sigma] = HIST_DIST(sample, v_ev, ...)
%
% Input:
%   sample         vector of positive, lognormally distributed values
%   v_ev           vector of values at which to evaluate PDF of `sample`
%   is_normalized  logical flag whether to normalize PDF of `sample`
%                  (default: true)
%
% Output:
%   pdf_max    position of maximum of probability density function (PDF)
%   err_max    absolute error of `pdf_max`
%   peak       value of lognormal PDF at `pdf_max`, normalized if
%              `is_normalized` is not false.
%   mu         first parameter of lognormal PDF
%   sigma      second parameter of lognormal PDF
%   err_mu     absolute error of `mu`
%   err_sigma  absolute error of `sigma`
%   pdf        PDF of `sample`, evaluated at values `v_ev`.
%              `pdf` is of the same size as `v_ev`.
%              If `v_ev` is empty or not given, `pdf` is an empty array.
%              If `is_normalized` is false, `pdf` is not normalized.
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

% Check input and output arguments
n_argin = nargin;
n_argout = abs(nargout);

if n_argin < 2
	is_pdf = false;
else
	is_pdf = true;
end

if n_argin < 3
	is_normalized = true;
else
	is_normalized = logical(is_normalized);
end

idx_max = 1;
idx_err_max = 0;
idx_peak = 0;
idx_pdf = 0;
idx_mu = 0;
idx_err_mu = 0;
idx_sigma = 0;
idx_err_sigma = 0;

if n_argout < 2
elseif is_pdf
	idx_err_max = 2;
	idx_pdf = 3;
	switch n_argout
		case 2
			idx_pdf = 2;
			idx_err_max = 0;
		case 3
		case 4
			idx_peak = 4;
		case 5
			idx_mu = 4;
			idx_sigma = 5;
		case 6
			idx_peak = 4;
			idx_mu = 5;
			idx_sigma = 6;
		case 7
			idx_mu = 4;
			idx_err_mu = 5;
			idx_sigma = 6;
			idx_err_sigma = 7;
		otherwise
			idx_peak = 4;
			idx_mu = 5;
			idx_err_mu = 6;
			idx_sigma = 7;
			idx_err_sigma = 8;
	end
else
	idx_err_max = 2;
	switch n_argout
		case 2
		case 3
			idx_peak = 3;
		case 4
			idx_mu = 3;
			idx_sigma = 4;
		case 5
			idx_peak = 3;
			idx_mu = 4;
			idx_sigma = 5;
		case 6
			idx_mu = 3;
			idx_err_mu = 4;
			idx_sigma = 5;
			idx_err_sigma = 6;
		otherwise
			idx_peak = 3;
			idx_mu = 4;
			idx_err_mu = 5;
			idx_sigma = 6;
			idx_err_sigma = 7;
	end
end

% Prepare and count sample
sample = sample(isfinite(sample));
n_sample = length(sample);

% Calculate sample mean and variance
sample_E = mean(sample);
sample_V = var(sample);

% Calculate distribution parameters
sigma = sqrt(log(sample_V / sample_E^2 + 1));
mu = log(sample_E) - sigma^2 * .5;

% Write position of sample distribution maximum
varargout{idx_max} = exp(mu - sigma^2);

% Write sample distribution parameters
if idx_mu
	varargout{idx_mu} = mu;
end
if idx_sigma
	varargout{idx_sigma} = sigma;
end

% Calculate errors of distribution properties
if idx_err_max || idx_err_mu || idx_err_sigma

	% Calculate error of sample mean and variance, using:
	%    Relative error of variance = sqrt(2/(n-1))
	% from: Squires (1971): "Meßergebnisse und ihre Auswertung", p. 37
	% (ISBN 3110036320)
	err_E = sqrt(sum((sample - sample_E).^2) / (n_sample * (n_sample-1)));
	err_V = sqrt(2 * sum((sample - sample_E).^2)) / (n_sample - 1);

	% Calculate derivatives of distribution parameters by
	% sample mean and variance
	Dsigma_DE = @(E,V) V / (V + E^2)^1.5;
	Dsigma_DV = @(E,V) 1 / (2 * (V + E^2) * sqrt(log(V / E^2 + 1)));
	Dmu_DE = @(E,V) (1 + V / (V + E^2)) / E;
	Dmu_DV = @(E,V) -1 / (2 * (V + E^2));

	% Calculate errors of distribution parameters
	err_mu = sqrt( (Dmu_DE(sample_E, sample_V) * err_E)^2 + ...
		(Dmu_DV(sample_E, sample_V) * err_V)^2 );
	err_sigma = sqrt( (Dsigma_DE(sample_E, sample_V) * err_E)^2 + ...
		(Dsigma_DV(sample_E, sample_V) * err_V)^2 );

	% Calculate and write error of distribution maximum
	if idx_err_max
		Dmax_Dmu = @(mu,sigma) exp(mu - sigma^2);
		Dmax_Dsigma = @(mu,sigma) -2 * sigma * exp(mu - sigma^2);

		varargout{idx_err_max} = sqrt( (Dmax_Dmu(mu, sigma) * err_mu)^2 + ...
			(Dmax_Dsigma(mu, sigma) * err_sigma)^2 );
	end

	% Write errors of distribution parameters
	if idx_err_mu
		varargout{idx_err_mu} = err_mu;
	end
	if idx_err_sigma
		varargout{idx_err_sigma} = err_sigma;
	end
end

% Calculate maximum of probability distribution
if idx_peak
	varargout{idx_peak} = pdf('lognormal', varargout{idx_max}, mu, sigma);
	if ~is_normalized
		varargout{idx_peak} = varargout{idx_peak} * n_sample;
	end
end

% Calculate probability distribution
if idx_pdf
	if isempty(v_ev) || ~is_pdf
		varargout{idx_pdf} = [];
	else
		varargout{idx_pdf} = pdf('lognormal', v_ev, mu, sigma);

		% Un-normalize probability distribution
		if ~is_normalized
			varargout{idx_pdf} = varargout{idx_pdf} * n_sample;
		end
	end
end
