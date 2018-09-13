% This script plots for the given conditions the event time histograms
% and the event time correlations.
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
cond_list = string({'ctrl','fas100','fas200','fas500','NP0.1','NP1','NP10','NP25','NP100','sts'});

for s = cond_list
	restrict_to_conditions = s;
% 	out_token = [s{:} '_'];
	plot_eventTimeHistograms
	plot_eventTimeCorrelations
	close all
end
