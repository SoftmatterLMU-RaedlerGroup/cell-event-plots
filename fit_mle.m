function [paramsOut, logL, params0, paramsFit, fmincon_output] = ...
	fit_mle(fun, indep, dep, varargin)
%fit_mle performs a maximum-likelihood estimation
%
% Input arguments:
%	fun			function handle, or cell array of function handles, of one
%				or more functions to fit
%	indep		independent variable ("x-values") to be used for fitting,
%				given as numeric array or cell array of numeric arrays. If
%				given as cell array, specify one or as many elements as fit
%				functions.
%	dep			dependent variable (measured values to which to fit `fun`).
%				For a single fit function, specify `dep` as numeric array.
%				Else, specify `dep` as cell array with as many numeric
%				elements as fit functions given.
%	varargin	"Name, Value" pairs (see below)
%
% Return values:
%	paramsOut	the parameter set that yields the largest likelihood
%	logL		vector of negative log-likelihood values in ascending order
%				(first value corresponds to `paramsOut`)
%	params0		the set of initial parameters used, sorted by descending
%				likelihood (first row corresponds to `paramsOut`)
%	paramsFit	the parameter sets yielded by the iterations, sorted by
%				descending likelihood (first row equals `paramsOut`)
%	fmincon_output	structure of all outputs from all calls to `fmincon`.
%				The fields are called like the corresponding output of
%				`fmincon` and contain vectors of the outputs, sorted by
%				descending likelihood (same ordering as for `logL`).
%				The field `exitflag` holds a sorted numerical vector, the
%				fields `output`, `lambda`, `grad` and `hessian` hold cell
%				vectors. The contents of the vector components is directly
%				taken from `fmincon`'s output.
%
% The function handles specified in `fun` must be of the form:
%	fun{i}(indep, parameters)
% where `indep` is the independent variable as given in `indep`, and
% `parameters` are the parameters taken by the function, either given as
% numerical vector or as comma-separated list.
%
% It is possible to specify multiple functions that use overlapping sets of
% parameters. To specify multiple functions
%	* give `fun` as a cell array of function handles,
%	* give `dep` as a cell array of as many value sets as functions,
%	* optionally give `indep` as a cell array of as many independent
%	  variables as functions,
%	* define a number of parameters via the keyword "Parameters", and
%	* define which function uses which parameters via the keyword
%	  "ParameterList".
%
% Possible "Name, Value" pairs:
%	ParameterBounds	Information on the parameters
%					Each column stands for one parameter. The rows have the
%					following order:
%					(1)	one row lower parameter bounds; default: -Inf
%					(2)	one row upper parameter bounds; default: +Inf
%					(3)	one row lower initial value bounds; default: -1
%					(4) one row upper initial value bounds; default: +1
%
%					Infinite values in (1) or (2) stand for unbound
%					parameters.
%					(1) and (2) must always be specified together. If
%					complete rows (1) and (2) have default values, they can
%					be replaced by one row of NaN.
%					(3) and (4) must always be specified together.
%					If (3) and (4) are specified, (1) and (2) must also be
%					specified.
%					A value in (1) resp. (3) must not be larger than the
%					corresponding value in (2) resp. (4).
%					(3) and (4) are only needed for generating initial
%					parameter sets via latin hypercube sampling. They are
%					adjusted to be consistent if the defaults are not
%					reasonable.
%
%	InitialParameters	Array of initial parameter values. The number of
%					columns gives the number of parameters, and the number
%					of rows determines the number of iterations.
%
%	ParameterList	Cell array of lists of indices for each function. This
%					must be specified when using multiple functions.
%					When using a single function, the parameter list can
%					be given as a numeric vector.
%
%	ParameterLog	Logical vector whose length equals the number of
%					parameters given. `true` means that the corresponding
%					parameter is fitted logarithmically; `false` means that
%					the parameter is fitted linearly.
%					When the same setting applies to all parameters, the
%					setting can be given as a scalar.
%					The default value is `false`.
%
%					For parameters to be treated logarithmically, the
%					values given by "ParameterBounds" or
%					"InitialParameters" are also interpreted as logarithmic
%					values. The return values are always the true values
%					and not their logarithms.
%
%	Iterations		The number of multistarts to be used. A latin hypercube
%					sample will be generated to ensure different initial
%					parameter values for each multistart iteration.
%					The default number of iterations is 10.
%					When an initial parameter value set is specified via
%					the "ParameterBounds" keyword, the number of iterations
%					is one if the "Iterations" keyword is not specified.
%
%	ObjectiveFcn	Function handle to a objective function that is used to
%					calculate the likelihood. The objective function values
%					of multiple functions are summed up.
%					The default objective function is the logarithm of a
%					normal distribution of the residuals.
%					A custom objective function must have the signature:
%						function(fun, indep, dep, parameters, doLogFit)
%					where `fun` is a function handle of the function to fit
%					which accepts the independent variable `indep` as only
%					argument, `dep` is the corresponding dependent
%					variable, `parameters` are the (non-logarithmic)
%					parameters used for this function evaluation, and
%					`doLogFit` is a scalar logical indicating whether to
%					perform a linear (`false`) or a logarithmic (`true`)
%					fit.
%					To use different objective functions for different fit
%					functions, specify the objective functions as cell
%					array of function handles. For empty cells, the default
%					objective function is used.
%
%	LogarithmicFit	Logical vector whose length is one or equals the number
%					of functions given. `true` means that the logarithm of
%					the corresponding function values are fitted to the
%					logarithm of the corresponding `dep` values; `false`
%					means that the (non-logarithmic) function values are
%					fitted to the (non-logarithmic) `dep` values.
%					If a scalar logical is given, the setting is used for
%					all functions.
%
%	Sigma			The standard deviation to be used by the default objective
%					function, given as double vector of length one or
%					equal to the number of functions. A scalar value is
%					applied to all functions. The default standard
%					deviation is 1.
%					An integer vector is interpreted as a vector of indices
%					to the parameters. Use this to fit the standard
%					deviation.
%					This value is only used in the default objective
%					function and ignored otherwise.
%
%	FminconOpt		Custom options to be used for `fmincon`. Specify the
%					options as cell array of Name-Value pairs as expected
%					by `optimoptions`, without the leading "fmincon"
%					argument.
%
% If not stated otherwise, any "Name, Value" pair is optional and can be
% omitted, in which case reasonable default values will be used.
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
%%

%% Read input
% Ensure that `fun` is a cell array of function handles
if isa(fun, 'function_handle')
	fun = {fun};
elseif iscell(fun)
	for iF = 1:numel(fun)
		if ~isa(fun{iF}, 'function_handle')
			error('Invalid function given at index %d', iF);
		end
	end
else
	error('No function given.');
end
nFuns = numel(fun);

% Ensure that `indep` is an array or a cell array
if isnumeric(indep)
	indep = {indep};
elseif iscell(indep)
	if numel(indep) ~= 1 && numel(indep) ~= nFuns
		error('Number of independent value sets does not match number of functions.');
	end
else
	error('Type of independent variable is not valid.');
end

% Ensure that `dep` is an array or a cell array
if isnumeric(dep)
	if nFuns > 1
		error(['When giving multiple functions, ' ...
			'specify as many dependent variable values as functions.']);
	end
	dep = {dep};
elseif iscell(dep)
	if numel(dep) ~= nFuns
		error('Number of dependent value sets does not match number of functions.');
	end
else
	error('Type of dependent variable is not valid.');
end


% Read out values from `varargin`
if mod(numel(varargin), 2)
	% `varargin` has odd number of elements
	error('`varargin` must be structured as: {key1, value1, key2, value2, ...}');
end

% Define default values
nParams = 0;
params = [];
params0 = [];
paramsList = [];
paramsLog = [];
sigma = [];
objective = [];
nIter = [];
fitLog = [];
fmincon_options = { ...
	'Algorithm', 'interior-point'; ...
	'MaxIter', 10000; ...
	'MaxFunEvals', 10000; ...
	'Display', 'notify' ...
	}';

% Iterate over `varargin`
iIn = 1;
while iIn < numel(varargin)

	if strcmpi(varargin{iIn}, 'ParameterBounds')
		% Find parameter bounds
		params = varargin{iIn+1};

	elseif strcmpi(varargin{iIn}, 'InitialParameters')
		% Get initial parameters
		params0 = varargin{iIn+1};

	elseif strcmpi(varargin{iIn}, 'ParameterList')
		% Get list of parameters to use
		paramsList = varargin{iIn+1};

	elseif strcmpi(varargin{iIn}, 'ParameterLog')
		% Find whether to fit parameters linearly or logarithmically
		paramsLog = varargin{iIn+1};

	elseif strcmpi(varargin{iIn}, 'Iterations')
		% Get number of iterations (fit rounds)
		nIter = varargin{iIn+1};

	elseif strcmpi(varargin{iIn}, 'ObjectiveFcn')
		% Get custom objective functions
		objective = varargin{iIn+1};

	elseif strcmpi(varargin{iIn}, 'LogarithmicFit')
		% Test whether to fit logarithmically
		fitLog = varargin{iIn+1};

	elseif strcmpi(varargin{iIn}, 'Sigma')
		% Find variance for default objective function
		sigma = varargin{iIn+1};

	elseif strcmpi(varargin{iIn}, 'FminconOpt')
		% Get custom options for `fmincon`
		fmincon_options = varargin{iIn+1};

	else
		% Throw error on unknown option
		error('Undefined option: %s', varargin{iIn});
	end

	% Increment `iIn`
	iIn = iIn + 2;
end

% Get parameters
% Format of parameter array:
%	* One column per parameter
%	* If parameter array is empty, use unbound parameters in [-1,+1]
%	* The next two rows are lower and upper bounds of parameters (default
%		value: ]-Inf,+Inf[; can be reduced to one row of NaN)
%	* The next two rows are lower and upper bounds for initial values of
%		unbound parameters (default: [-1,+1]; ignored for finite bounds)

% Test if initial parameters are given
if ~isempty(params0)
	nParams = size(params0, 2);
	if any(isnan(params0(:)))
		error('NaN given as initial parameter value.')
	end
end

% Get number of parameters
if nParams == 0
	if isempty(params)
		% Test if parameter array is empty
		for iF = 1:nFuns
			nParams = max(nParams, nargin(fun{iF}) - 1);
		end
		params = NaN(1, nParams);
	else
		nParams = size(params, 2);
	end
elseif isempty(params)
	params = NaN(1, nParams);
elseif size(params, 2) ~= nParams
	error('Number of initial parameter values does not match number of bounds.')
end

% Extract lower/upper bound vector (i.e. all NaN)
paramLowerBnd = -Inf(1, nParams);
paramUpperBnd = Inf(1, nParams);
offsetPar = 0;

if (size(params,1) >= offsetPar+1) && any(isnan(params(offsetPar+1,:)))
	if all(isnan(params(offsetPar+1,:)))
		offsetPar = offsetPar + 1;
	else
		error('Parameter bounds must be either all NaN (unbound) or no NaN.');
	end

elseif size(params, 1) >= offsetPar + 2
	paramLowerBnd = params(offsetPar + 1, :);
	paramUpperBnd = params(offsetPar + 2, :);
	offsetPar = offsetPar + 2;

	if any(isnan(paramUpperBnd))
		error('Parameter bounds must be either all NaN (unbound) or no NaN.');
	elseif any(paramLowerBnd > paramUpperBnd)
		error('Lower parameter bounds must not be larger than upper bounds.');
	end
end

% Extract initial parameter value bounds
paramLower0 = paramLowerBnd;
paramUpper0 = paramUpperBnd;

if (size(params,1) >= offsetPar+2)
	idxLo = isinf(paramLower0);
	paramLower0(idxLo) = params(offsetPar + 1, idxLo);
	idxUp = isinf(paramUpper0);
	paramUpper0(idxUp) = params(offsetPar + 2, idxUp);

	if any(~isfinite([paramLower0,paramUpper0]))
		error('Initial parameter value bounds must be finite.');
	elseif any(paramLower0 > paramUpper0)
		error('Lower initial parameter bounds must not be larger than upper ones.');
	elseif any(paramLower0 < paramLowerBnd | paramUpper0 > paramUpperBnd)
		error('Initial parameter bounds must not exceed general parameter bounds.');
	end

else
	idxLo = isinf(paramLower0);
	idxUp = isinf(paramUpper0);

	paramLower0(idxLo) = -1;
	paramUpper0(idxUp) = 1;

	idx_temp = idxLo & ~idxUp & paramUpper0 < 0;
	paramLower0(idx_temp) = paramUpper0(idx_temp) - 1;

	idx_temp = idxUp & ~idxLo & paramLower0 > 0;
	paramUpper0(idx_temp) = paramLower0(idx_temp) + 1;
end

% Test list of parameters
paramsIdx = {};
if isempty(paramsList)
	for iF = 1:nFuns
		nIn = nargin(fun{iF}) - 1;
		if nIn < 1
			error('Too few input arguments in function %d', iF);
		elseif nIn > nParams
			error('Parameter index (%d) exceeds parameter number.', nIn);
		end
		paramsIdx{iF} = 1:nIn;
	end
elseif iscell(paramsList)
	if numel(paramsList) ~= nFuns
		error('Number of parameter lists (%d) does not match number of functions (%d).', ...
			numel(paramsList), nFuns);
	end

	for iF = 1:nFuns
		nIn = nargin(fun{iF}) - 1;
		nInGiven = numel(paramsList{iF});
		if nInGiven < nIn
			error('Too few parameters given for function %d: %d required, %d given.', ...
				iF, nIn, nInGiven);
		elseif nInGiven > nIn && nIn ~= 1
			error('Too many input parameters given for function %d', iF);
		end
		paramsIdx{iF} = paramsList{iF};
	end
elseif isnumeric(paramsList)
	if nFuns > 1
		error('Parameter indices for more than one function must be given as cell array.');
	end
	paramsIdx = {paramsList(:)};
else
	error('List of parameter indices has invalid format.');
end

% Test if parameters are linear or logarithmical
if ~isempty(paramsLog)
	if ~isvector(paramsLog)
		paramsLog = paramsLog(:);
	end

	if ~islogical(paramsLog)
		error('The ParameterLog option requires a vector suitable for indexing the parameters.');
	elseif length(paramsLog) == 1
		paramsLog(1:nParams) = paramsLog;
	elseif length(paramsLog) ~= nParams
		error('The ParameterLog vector length does not match the parameter number.');
	end
else
	paramsLog = false(1, nParams);
end

% Get number of iterations and initial parameter sets
if isempty(nIter)
	if isempty(params0)
		nIter = 10;
	else
		nIter = size(params0, 1);
	end
elseif ~isnumeric(nIter) || nIter < 1
	error('The Iterations option requires a positive integer number.')
elseif nIter < size(params0, 1)
	nIter = size(params0, 1);
	warning('Setting number of iterations to number of initial parameter sets.');
end

if nIter - size(params0, 1) > 0
	% Create normalized latin hypercube sample
	lhs = lhsdesign(nIter - size(params0, 1), nParams, 'criterion', 'none');
	lhs = paramLower0 + (paramUpper0 - paramLower0) .* lhs;
	params0 = [ params0; lhs ];
end

% Test whether to fit logarithmically
if isempty(fitLog)
	fitLog = false(1, nFuns);
else
	if ~islogical(fitLog)
		error('Option "LogarithmicFit" requires a logical vector.');
	end
	if ~isvector(fitLog)
		fitLog = fitLog(:)';
	end
	if length(fitLog) ~= 1 && length(fitLog) ~= nFuns
		error(['Option "LogarithmicFit" requires a logical scalar or ' ...
			'vector whose length matches the number of functions.']);
	elseif length(fitLog) == 1
		fitLog(1:nFuns) = fitLog;
	end
end

% Test value for sigma
if ~isempty(sigma)

	if isa(sigma, 'double')
		if sigma <= 0
			error('The standard deviation must not be equal or less than zero.');
		end
	elseif isinteger(sigma)
		if any(sigma < 1 | sigma > nParams)
			error('The standard deviation reference must be a parameter index.');
		end
	else
		error(['The standard deviation of the default objective ' ...
			'function must be of type integer or double.']);
	end

	if ~isvector(sigma)
		sigma = sigma(:).';
	end

	if length(sigma) ~= 1 && length(sigma) ~= nFuns
		error(['The length of the standard deviation of the default objective ' ...
			'function must be 1 or the length of the numer of functions.']);
	elseif length(sigma) == 1
		sigma(2:nFuns) = sigma;
	end
else
	sigma = ones(1, nFuns);
end

% Test presence of objective functions
if ~isempty(objective)
	if isa(objective, 'function_handle')
		objective = repmat({objective}, 1, nFuns);
	elseif iscell(objective)
		if numel(objective) ~= nFuns && numel(objective) ~= 1
			error(['Number of objective functions must be one or match '...
				number of fit functions.']);
		end
		for iF = 1:numel(objective)
			if ~isa(objective{iF}, 'function_handle') ...
					&& ~isempty(objective{iF})
				error('Objective functions must be function handles.');
			end
		end
	else
		error(['Option "ObjectiveFcn" requires a function handle or a ' ...
			'cell array of function handles.']);
	end
else
	objective = cell(1, nFuns);
end

%% Perform multistart fit
logL = NaN(nIter, 1);
paramsFit = NaN(nIter, nParams);
exitflags = NaN(nIter, 1);
outputs = cell(nIter, 1);
lambdas = cell(nIter, 1);
grads = cell(nIter, 1);
hessians = cell(nIter, 1);
fmincon_options = optimoptions('fmincon', fmincon_options{:});

for iIter = 1:nIter
	[paramsFit(iIter,:), logL(iIter), exitflags(iIter), outputs{iIter}, ...
		lambdas{iIter}, grads{iIter}, hessians{iIter}] = ...
		fmincon(@applyObjective, params0(iIter,:), [], [], [], [], ...
		paramLowerBnd, paramUpperBnd, [], fmincon_options);
end

%% Prepare output
[logL, logL_idx] = sort(logL);

paramsFit = paramsFit(logL_idx,:);
paramsFit(:,paramsLog) = 10.^paramsFit(:,paramsLog);

paramsOut = paramsFit(1,:);

params0 = params0(logL_idx,:);
params0(:,paramsLog) = 10.^params0(:,paramsLog);

fmincon_output = struct(...
	'exitflag', exitflags(logL_idx), ...
	'output', outputs(logL_idx), ...
	'lambda', lambdas(logL_idx), ...
	'grad', grads(logL_idx), ...
	'hessian', hessians(logL_idx));


%% Function called by `fmincon`
function logL = applyObjective(params)
%applyObjective is passed to `fmincon` and determines the quality of fit

	% Initialize log-likelihood
	logL = 0;

	% Get indices for independent variable
	if numel(indep) == 1
		iIndep = ones(1, nFuns);
	else
		iIndep = 1:nFuns;
	end

	% Convert parameters from log-space to linear space
	params(paramsLog) = 10.^params(paramsLog);

	% Get standard deviation for default objective function
	if isinteger(sigma)
		theseSigma = params(sigma);
	else
		theseSigma = sigma;
	end

	% Iterate over functions
	for jF = 1:nFuns

		% Get the parameters for this function
		thisParams = params(paramsIdx{jF});
		if nargin(fun{jF}) - 1 == 1 && length(thisParams) > 1
			thisParams = {thisParams};
		else
			thisParams = num2cell(thisParams);
		end

		% Build a function handle for the given parameter set
		thisFun = @(independent) fun{jF}(independent, thisParams{:});

		% Calculate negative log-likelihood
		if isempty(objective{jF})
			% No objective function specified; use default
			logL = logL + defaultObjective(thisFun, indep{iIndep(jF)}, ...
				dep{jF}, fitLog(jF), theseSigma(jF));
		else
			% Use objective function specified by user
			logL = logL + objective{jF}(thisFun, indep{iIndep(jF)}, ...
				dep{jF}, params, fitLog(jF));
		end
	end

end % end of applyObjective

end % end of `fit_mle`

%% Auxiliary functions

function logL = defaultObjective(thisFun, thisIndep, thisDep, doLogFit, sigma)
%defaultObjective is the default objective function that returns a
%log-likelihood

	% Evaluate the function
	thisVals = thisFun(thisIndep);

	% Calculate the residuals
	if doLogFit
		residuals = log(thisDep) - log(thisVals);
	else
		residuals = thisDep - thisVals;
	end
	residuals = residuals(:);

	% Calculate the log-likelihood (more specifically: -2*log(L)
	QstdDev = sigma^2;
	logL = length(residuals) * log(QstdDev * 2 * pi) + ...
		sum(residuals.^2) / QstdDev;
end % end of `defaultObjective`
