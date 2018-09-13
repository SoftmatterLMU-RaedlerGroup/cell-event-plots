%test_ms is a mean shift test
% Required: d (dim-by-occurrence matrix)
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

% Prepare values
S = d;
h = 3;
doLoop = true;
iRun = 0;

% Plot initial sample
ax = plotKDE(S, h);
ax = plotS(S, S, h, ax);
title(ax, 'Initial sample')

% Perform mean shift
Told = S;
while doLoop
	% Calculate mean shift
	T = kde_mean(Told, S, h);

	% Test for finish
	if all( sqrt(sum((T - Told).^2, 1)) < 0.02 )
		% Iteration finished
		doLoop = false;
	else
		Told = T;
	end

	% Plot
	iRun = iRun + 1;
	if ~doLoop %|| mod(iRun, 20) == 0
		ax = plotKDE(S, h);
		ax = plotS(T, S, h, ax);
		title(ax, sprintf('Run %d', iRun))
	end
end

% Find clusters
labels = zeros(1, size(T, 2));
new_lbl = 1;
doLoop = true;
while doLoop
	% Find data point without label
	iFree = find(~labels, 1);
	if isempty(iFree)
		doLoop = false;
		break
	end

	% Find similar data points
	iNew = sqrt(sum((T - T(:,iFree)).^2, 1)) < 1/6;
	lbls = unique(labels(iNew));
	lbls(lbls == 0) = [];

	if isempty(lbls)
		% No label assigned yet, assign new label
		labels(iNew) = new_lbl;
		new_lbl = new_lbl + 1;
	elseif length(lbls) == 1
		% Other label found, use it
		labels(iNew) = lbls;
	else
		% More than one other label found, use smallest
		nl = min(lbls);
		labels(iNew) = nl;

		% Assing new label to all connected clusters
		for l = lbls(lbls ~= nl)
			labels(labels == l) = nl;
		end
	end
end

% Reduce cluster number
lbl_list = unique(labels);
new_labels = zeros(size(labels));
nLblOccur = zeros(size(lbl_list));
nLbls = length(lbl_list);
for iL = 1:nLbls
	thisIdx = ismember(labels, lbl_list(iL));
	new_labels(thisIdx) = iL;
	nLblOccur(iL) = sum(thisIdx);
end
labels = new_labels;

% Build colormap for clusters
clr = zeros(nLbls, 3);
nLblMoreThanOnce = sum(nLblOccur > 1);
cm = lines(nLblMoreThanOnce);
iLblMoreThanOnce = 1;
for iL = 1:nLbls
	if nLblOccur(iL) > 1
		clr(iL, :) = cm(iLblMoreThanOnce, :);
		iLblMoreThanOnce = iLblMoreThanOnce + 1;
	end
end

% Plot colored clusters
f = figure;
ax = axes(f);
hold(ax, 'on');
line(ax, [0, 30], [0, 30], 'color', 'k');
for iS = 1:size(S, 2)
	plot(ax, S(1,iS), S(2, iS), '.', 'color', clr(labels(iS), :))
end
ax.XLim = [0,30];
ax.YLim = [0,30];
xlabel(ax, 't1 [h]')
ylabel(ax, 't2 [h]')

% Plot with profile
ax = plotKDE(S, h);
for iS = 1:size(S, 2)
	plot3(ax, S(1,iS), S(2,iS), kde(S(:,iS), S, h), '.', ...
		'color', clr(labels(iS), :))
end
xlabel(ax, 't1 [h]')
ylabel(ax, 't2 [h]')


function ax = plotS(T, S, h, ax)
%plotS plots the sample points of a KDE
p = zeros(1, size(T, 2));
for iT = 1:size(T, 2)
	p(iT) = kde(T(:,iT), S, h);
end
plot3(ax, T(1,:) , T(2,:), p, '.r')
end


function ax = plotKDE(S, h)
%plotKDE plots a KDE surface
	% Make density profile
	x = linspace(0, 30, 100);
	y = linspace(0, 30, 100);
	z = zeros(length(y), length(x));

	lbl1 = 't1 [h]';
	lbl2 = 't2 [h]';

	for ix = 1:length(x)
		for iy = 1:length(y)
			z(iy,ix) = kde([x(ix);y(iy)], S, h);
		end
	end

	% Plot
	f = figure;
	ax = axes(f);
	hold(ax, 'on')
	surf(ax, x, y, z, 'EdgeColor', 'none', 'FaceColor', 'interp')
	xlabel(ax, lbl1)
	ylabel(ax, lbl2)
end


function k = kde(T, S, h)
%kde calculates a kernel density estimation of data S at points T with
% bandwidth h
	k = zeros(1, size(T, 2));
	for iT = 1:size(T, 2)
		k(iT) = sum(K(sqrt(sum((S - T(:,iT)).^2, 1)) / h));
	end
end


function k = kde_mean(T, S, h)
%kde calculates a sample mean of data S at points T with bandwidth h
	k = zeros(size(T));
	for iT = 1:size(T, 2)
		k(:,iT) = sum( K(sqrt(sum((S - T(:,iT)).^2, 1)) / h) .* S, 2 ) ...
			/ kde(T(:,iT), S, h);
	end
end


% function k = K(x)
% %K Gaussian kernel
% k = exp(-x.^2);
% end


function k = K(x)
%K epanechnikov kernel
	k = zeros(size(x));
	idx = abs(x) <= 1;
	k(idx) = 0.75 * (1 - x(idx).^2);
end