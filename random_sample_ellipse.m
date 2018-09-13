%random_sample_ellipse plots random data with an asymmetric ellipse.
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

%% Define paramters
nPts = 200;
muX = 1.5;
muY = 1.5;
sigmaX1 = .5;
sigmaX2 = 1.2;
sigmaY = .3;
theta = 50;
maxX = 4;
maxY = 4;

%% Generate random sample
rng('shuffle')
seed = randi(bitcmp(uint32(0)));
rng(seed)

x = [ -abs(normrnd(0, sigmaX1, 1, nPts/2)), abs(normrnd(0, sigmaX2, 1, nPts/2)) ];
y = normrnd(0, sigmaY, 1, nPts);

R_xy = [cosd(theta), -sind(theta); sind(theta), cosd(theta)];

P = R_xy * [x;y] + [muX;muY];

[ell, x0, y0, ~, ~, ~, ~, ~, info] = error_potato(P', 'large');
ell = [ell(:,1)+x0, ell(:,2)+y0];

R_ell = [cos(info.angle), -sin(info.angle); sin(info.angle), cos(info.angle)];
saXpos = R_ell * [0, info.stdXpos; 0, 0] + [x0; y0];
saXneg = R_ell * [0, -info.stdXneg; 0, 0] + [x0; y0];
saYpos = R_ell * [0, 0; 0, info.stdYpos] + [x0; y0];
saYneg = R_ell * [0, 0; 0, -info.stdYneg] + [x0; y0];

%% Plot
fh = figure;
ax = axes(fh);
hold(ax, 'on')

plot(ax, P(1,:), P(2,:), '.b')
plot(ax, saXpos(1,:), saXpos(2,:), '-r')
plot(ax, saXneg(1,:), saXneg(2,:), '-m')
plot(ax, saYpos(1,:), saYpos(2,:), '-g')
plot(ax, saYneg(1,:), saYneg(2,:), '-c')
plot(ax, ell(:,1), ell(:,2), '-k')
plot(ax, x0, y0, '+k')

xlabel(ax, 't_1')
ylabel(ax, 't_2')
title(sprintf('Random seed: %d', seed))

ax.XLim = [0, maxX];
ax.YLim = [0, maxY];
pbaspect(ax, [1 1 1]);
ax.XTick = [];
ax.YTick = [];
ax.Box = 'on';

%% Save
outfilename = fullfile('out', [getTime, '_random_potato.pdf']);

fh.PaperUnits = 'centimeters';
fh.PaperSize = [10 10];
fh.PaperPosition = [0 0 fh.PaperSize];

print(fh, '-painters', outfilename, '-dpdf')