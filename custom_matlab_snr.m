function [r, noisePow] = custom_matlab_snr(varargin)
%SNR    Signal to Noise Ratio
%   R = SNR(X, Y) computes the signal to noise ratio (SNR) in dB, by
%   computing the ratio of the summed squared magnitude of the signal, X,
%   to the summed squared magnitude of the noise, Y, where Y has the same
%   dimensions as X.  Use this form of SNR when your input signal is not
%   sinusoidal and you have an estimate of the noise.
%
%   R = SNR(X) computes the signal to noise ratio (SNR), in dBc, of
%   the real sinusoidal input signal X.  The computation is performed over
%   a periodogram of the same length as the input using a Kaiser window
%   and excludes the power of first six harmonics (including the
%   fundamental).
%
%   R = SNR(X, Fs, N, tf, max_f_dc) computes the signal to noise ratio (SNR) in 
%   dBc, of
%   the real sinusoidal input signal, X, with sampling rate, Fs, and number
%   of harmonics, N, to exclude from computation when computing SNR.  The
%   default value of Fs is 1.  The default value of N is 6 and includes the
%   fundamental frequency. If tf is porvided, use it as fundamental frequency.
%   If tf's length is large than 1, the first frequency will be seen as target 
%   frequency, and the others will be seen as harmonics.
%   If max_f_dc, the maximum frequency of DC component, 
%   is provided, the spectrum in frequency interval [0, max_f_dc] will be
%   replaced with the median spectrum before counting noise power.
%
%   R = SNR(Pxx, F, 'psd', tf, max_f_dc) specifies the input as a one-sided 
%   PSD estimate,
%   Pxx, of a real signal.   F is a vector of frequencies that corresponds
%   to the vector of Pxx estimates.  The computation of noise excludes the
%   first six harmonics (including the fundamental).  If tf is porvided, use 
%   it as fundamental frequency. If tf's length is large than 1, the first
%   frequency will be seen as target frequency, and the others will be seen
%   as harmonics.
%   If max_f_dc, the maximum frequency of DC component, 
%   is provided, the spectrum in frequency interval [0, max_f_dc] will be
%   replaced with the median spectrum before counting noise power.
%
%   R = SNR(Pxx, F, N, 'psd', tf, max_f_dc) specifies the number of harmonics, 
%   N, to
%   exclude when computing SNR.  The default value of N is 6 and includes
%   the fundamental frequency.  If tf is porvided, use it as fundamental 
%   frequency. If tf's length is large than 1, the first frequency will be 
%   seen as target frequency, and the others will be seen as harmonics.
%   If max_f_dc, the maximum frequency of DC component, 
%   is provided, the spectrum in frequency interval [0, max_f_dc] will be
%   replaced with the median spectrum before counting noise power.
%
%   R = SNR(Sxx, F, RBW, 'power', tf, max_f_dc) specifies the input as a 
%   one-sided power
%   spectrum, Sxx, of a real signal.  RBW is the resolution bandwidth over
%   which each power estimate is integrated. If tf is porvided, use 
%   it as fundamental frequency. If tf's length is large than 1, the first
%   frequency will be seen as target frequency, and the others will be seen
%   as harmonics. If max_f_dc, the maximum frequency of DC component, 
%   is provided, the spectrum in frequency interval [0, max_f_dc] will be
%   replaced with the median spectrum before counting noise power.
%
%   R = SNR(Sxx, F, RBW, N, 'power', tf, max_f_dc) specifies the number of 
%   harmonics, N,
%   to exclude when computing SNR.  The default value of N is 6 and
%   includes the fundamental frequency. If tf is porvided, use 
%   it as fundamental frequency. If tf's length is large than 1, the first
%   frequency will be seen as target frequency, and the others will be seen
%   as harmonics. If max_f_dc, the maximum frequency of DC component, 
%   is provided, the spectrum in frequency interval [0, max_f_dc] will be
%   replaced with the median spectrum before counting noise power.
%
%   [R, NOISEPOW] = SNR(...) also returns the total noise power of the non-
%   harmonic components of the signal.
%
%   [...] = SNR(...,'aliased') also excludes harmonics of the fundamental
%   aliased into the Nyquist range.  Use this option when the input signal
%   is undersampled.  If unspecified or 'omitaliases' is used instead,
%   harmonics of the fundamental frequency beyond Nyquist are treated as
%   noise.  This works for all signatures listed above except SNR(X, Y).
%
%   SNR(...) with no output arguments plots the spectrum of the signal and
%   annotates the fundamental, DC component, harmonics and noise.  The DC
%   component is removed before computing SNR.  This works for all
%   signatures listed above except SNR(X, Y).
%
%   % Example 1:
%   %   Compute the SNR of a 2 second 20ms rectangular pulse sampled at
%   %   10 kHz in the presence of gaussian noise
%
%   Tpulse = 20e-3; Fs = 10e3;
%   x = rectpuls((-1:1/Fs:1),Tpulse);
%   y = 0.00001*randn(size(x));
%   s = x + y;
%   pulseSNR = snr(x,s-x)
%
%   % Example 2:
%   %   Plot the SNR of a 2.5 kHz distorted sinusoid sampled at 48 kHz
%   load('sineex.mat','x','Fs');
%   snr(x,Fs)
%
%   % Example 3:
%   %   Generate the periodogram of a 2.5 kHz distorted sinusoid sampled
%   %   at 48 kHz and measure the SNR (in dB)
%   load('sineex.mat','x','Fs');
%   w = kaiser(numel(x),38);
%   [Sxx, F] = periodogram(x,w,numel(x),Fs,'power');
%
%   % Measure SNR on the power spectrum
%   rbw = enbw(w,Fs);
%   sineSNR = snr(Sxx,F,rbw,'power')
%
%   % annotate the spectrum
%   snr(Sxx,F,rbw,'power')
%
%   See also SINAD THD SFDR TOI.

%   Copyright 2013-2019 The MathWorks, Inc.

%#codegen

narginchk(1, 7);   %%%%%%%%%%%%%%%%%%%

inputArgs = cell(1,length(varargin));
[inputArgs{:}] = convertStringsToChars(varargin{:});

% handle canonical definition of SNR
if nargin == 2 && isequal(size(inputArgs{1}), size(inputArgs{2}))
    [r, noisePow] = sampleSNR(inputArgs{1}, inputArgs{2});
    return
end

if coder.target('MATLAB') % for MATLAB execution
    
    % look for psd or power window compensation flags
    [estType, Args, ~] = signal.internal.psdesttype({'psd','power','time'},'time',inputArgs);
    
    % allow undersampled harmonics when 'aliased' is specified
    [harmType, Args] = getmutexclopt({'aliased','omitaliases'},'omitaliases',Args);
    
else % for code generation
    
    % look for psd or power window compensation flags and find their index
    % in the input argument list
    [estType, estIdx] = localInputParser({'psd','power','time'},'time',inputArgs);
    
    % look for aliased or omit aliases options and find their index in the
    % input argument list
    [harmType, harmIdx] = localInputParser({'aliased','omitaliases'},'omitaliases',inputArgs);
    
    % remove the above string options, if present, from the input argument
    % list
    if (estIdx~=0) && (harmIdx~=0)
        Args = cell(1,length(inputArgs)-2);
        if estIdx < harmIdx
            [Args{:}] = inputArgs{[1:(estIdx-1), (estIdx+1):(harmIdx-1), (harmIdx+1):end]};
        else
            [Args{:}] = inputArgs{[1:(harmIdx-1), (harmIdx+1):(estIdx-1), (estIdx+1):end]};
        end
    elseif (estIdx~=0) && (harmIdx==0)
        Args = cell(1,length(inputArgs)-1);
        [Args{:}] = inputArgs{[1:(estIdx-1), (estIdx+1):end]};
    elseif (estIdx==0) && (harmIdx~=0)
        Args = cell(1,length(inputArgs)-1);
        [Args{:}] = inputArgs{[1:(harmIdx-1), (harmIdx+1):end]};
    else
        Args = inputArgs;
    end
    
end

% check for unrecognized strings
chknostropts(Args{:});

% plot if no arguments are specified
plotType = distplottype(nargout, estType);

switch estType
    case 'psd'
        [r, noisePow] = psdSNR(plotType, harmType, Args{:});
    case 'power'
        [r, noisePow] = powerSNR(plotType, harmType, Args{:});
    case 'time'
        [r, noisePow] = timeSNR(plotType, harmType, Args{:});
end
end
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function [r, noisePow] = sampleSNR(x, y)

validateattributes(x,{'double','single'},{'real','finite','vector'}, ...
    'mysnr','X',1);

validateattributes(y,{'double','single'},{'real','finite','vector'}, ...
    'mysnr','Y',1);

signalPow = sum(x(:).^2);
noisePow  = sum(y(:).^2);
r = 10 * log10(signalPow / noisePow);
end
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function [r, noisePow] = timeSNR(plotType, harmType, x, fs, nHarm, tf, max_f_dc)

% force column vector before checking attributes
if length(x) == numel(x)
    colX = x(:);
else
    colX = x;
end

validateattributes(colX,{'double','single'},{'real','finite','vector'}, ...
    'mysnr','X',1);

if nargin > 3
    validateattributes(fs,{'numeric'},{'real','finite','scalar','positive'}, ...
        'mysnr','Fs',2);
    fsScalar = double(fs(1));
else
    fsScalar = 1;
end

if nargin > 4
    validateattributes(nHarm,{'numeric'},{'integer','finite','positive','scalar','>=',1}, ...
        'mysnr','N',3);  %%%%%%%%%%%%%%%%%%%%
    nHarmScalar = double(nHarm(1));
else
    nHarmScalar = 6;
end

%%%%%%%%%%%%%%%%%%%
if nargin < 6
    tf = [];
end

if nargin < 7
    max_f_dc = [];
end
%%%%%%%%%%%%%%%%%%%%

% remove DC component
colX = colX(:,1) - mean(colX(:,1));

n = length(colX);

% use Kaiser window to reduce effects of leakage
w = kaiser(n,38);
rbw = enbw(w,fsScalar);
[Pxx, F] = periodogram(colX,w,n,fsScalar,'psd');

% specify sample rate for aliased harmonics (assume same as Fs)
if strcmp(harmType, 'aliased')
    aliasFs = fsScalar;
else
    aliasFs = NaN;
end

[r, noisePow] = computeSNR(plotType, Pxx, F, rbw, nHarmScalar, aliasFs, tf, max_f_dc);  %%
end
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

function [r, noisePow] = powerSNR(plotType, harmType, Sxx, F, rbw, nHarm, tf, max_f_dc)  %%

validateattributes(Sxx, {'single','double'},...
    {'real','finite','nonnegative','vector'}, 'mysnr', 'Sxx', 1);
validateattributes(F, {'single','double'},...
    {'real','finite','vector','increasing','numel',numel(Sxx)}, ...
    'mysnr', 'F', 2);

coder.internal.assert(F(1)==0,'signal:snr:MustBeOneSidedSxx');

% ensure specified RBW is larger than a bin width
df = mean(diff(F));

validateattributes(rbw,{'numeric'},{'real','finite','positive','scalar','>=',df}, ...
    'mysnr','RBW',3);
rbwScalar = double(rbw(1));

if nargin > 5
    validateattributes(nHarm,{'numeric'},{'integer','finite','positive','scalar','>=',1}, ...
        'mysnr','N',4);  %%%%%%%%%%%%%%%%%%%%
    nHarmScalar = double(nHarm(1));
else
    nHarmScalar = 6;
end

%%%%%%%%%%%%%%%%%%%
if nargin < 7
    tf = [];
end

if nargin < 8
    max_f_dc = [];
end
%%%%%%%%%%%%%%%%%%%%

% assume frequency transform is from an even-length input when
% undersampling harmonics (i.e. F(end) = Fs/2)
if strcmp(harmType, 'aliased')
    aliasFs = 2*F(end);
else
    aliasFs = NaN;
end

[r, noisePow] = computeSNR(plotType, Sxx/rbwScalar, F, rbwScalar, nHarmScalar, aliasFs, tf, max_f_dc); %%
end
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function [r, noisePow] = psdSNR(plotType, harmType, Pxx, F, nHarm, tf, max_f_dc)  %%

validateattributes(Pxx, {'single','double'},...
    {'real','finite','nonnegative','vector'}, 'mysnr', 'Pxx', 1);
validateattributes(F, {'single','double'},...
    {'real','finite','vector','increasing','numel',numel(Pxx)}, ...
    'mysnr', 'F', 2);

coder.internal.assert(F(1)==0,'signal:snr:MustBeOneSidedPxx');

% use the average bin width
df = mean(diff(F));

if nargin > 4
    validateattributes(nHarm,{'numeric'},{'integer','finite','positive','scalar','>=',1}, ...
        'mysnr','N',4);  %%%%%%%%%%%
    nHarmScalar = double(nHarm(1));
else
    nHarmScalar = 6;
end

%%%%%%%%%%%%%%%%%%%
if nargin < 6
    tf = [];
end
if nargin < 7
    max_f_dc = [];
end
%%%%%%%%%%%%%%%%%%%%

% assume frequency transform is from an even-length input when
% undersampling harmonics (i.e. F(end) = Fs/2)
if strcmp(harmType, 'aliased')
    aliasFs = 2*F(end);
else
    aliasFs = NaN;
end

[r, noisePow] = computeSNR(plotType, Pxx, F, df, nHarmScalar, aliasFs, tf, max_f_dc);  %%
end
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function [r, noisePow] = computeSNR(plotType, Pxx, F, rbw, nHarm, aliasFs, tf, max_f_dc)  %%

%%%%%%%%%%%%%%%%%%%
nGivenHarm = 0;
if nargin < 7
    tf = [];
else
    nGivenHarm = numel(tf) - 1;
end
if nargin < 8
    max_f_dc = [];
end
%%%%%%%%%%%%%%%%%%%%

% save a copy of the original PSD estimates
origPxx = Pxx;

% pre-allocate harmonic table
psdHarmPow = NaN(nHarm+nGivenHarm,1);
psdHarmFreq = NaN(nHarm+nGivenHarm,1);
harmIdx = NaN(nHarm+nGivenHarm, 2);

% bump DC component by 3dB and remove it.
Pxx(1) = 2*Pxx(1);
Pxx_bump_dc = Pxx;  %%%

[~, ~, ~, iLeft, iRight] = signal.internal.getToneFromPSD(Pxx_bump_dc, F, rbw, 0);  %%%
if ~isempty(iLeft) && ~isempty(iRight)
    %%%%%%%
    if isempty(tf) || F(iRight) < tf(1)
        Pxx(iLeft(1):iRight(1)) = 0;    
    end 
    %%%%%%%
end
%%%%%%%%
if ~isempty(tf) && F(iRight) >= tf(1)  
    % target frequency has been recognized as DC
    Pxx(1) = 0;
    dcIdx = [1; 1];
else
    dcIdx = [iLeft; iRight];
end % if ~isempty(tf) && F(iRight) >= tf(1)
%%%%%%%%


if isempty(tf)  %%%%%%%%%%%%%%
    % get an estimate of the actual frequency / amplitude, then remove it.
    [Pfund, Ffund, iFund, iLeft, iRight] = signal.internal.getToneFromPSD(Pxx_bump_dc, F, rbw);  %%%
else   %%%%%%%%%%%%%%%%%%%%
    [Pfund, Ffund, iFund, iLeft, iRight] = signal.internal.getToneFromPSD(Pxx_bump_dc, F, rbw, tf(1)); %%
end   %%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%
%  force remove dc
if ~isempty(max_f_dc)
    max_f_dc = min(F(iLeft), max_f_dc);
    max_f_dc = max(max_f_dc, F(dcIdx(2)));
    for fi = 1:length(F)
        if F(fi) < max_f_dc
            Pxx(fi) = 0;
        else  % out of dc frequency range
            dcIdx(2) = fi - 1;
            break;
        end % of if F(fi) < max_f_dc
    end % of for fi = 1:length(F)

end % if ~isempty(max_f_dc)

%%%%%%%%%%%%%%%%

[psdHarmPow(1), psdHarmFreq(1)] = idx2psddb(Pxx, F, iFund);
if ~isempty(iLeft) && ~isempty(iRight)
    Pxx(iLeft(1):iRight(1)) = 0;
    harmIdx(1, :) = [iLeft(1); iRight(1)];
end

% remove harmonic content
for i=2:nHarm
    toneFreq = i*Ffund;
    if isfinite(aliasFs)
        toneFreq = aliasToNyquist(i*Ffund, aliasFs);
    end
    [harmPow, ~, iHarm, iLeft, iRight] = signal.internal.getToneFromPSD(Pxx_bump_dc, F, rbw, toneFreq);  %%%
    [psdHarmPow(i), psdHarmFreq(i)] = idx2psddb(Pxx, F, iHarm);
    % obtain local maximum value in neighborhood of bin
    if ~isnan(harmPow)
        if ~isempty(iLeft) && ~isempty(iRight)
            % remove the power of this tone
            Pxx(iLeft(1):iRight(1)) = 0;
            harmIdx(i, :) = [iLeft(1); iRight(1)];
        end
    end
end

%%%%%%%%%%%%%%%%
% remove given harmonic content
if length(tf) > 1
    for i = 2:numel(tf)
        toneFreq = tf(i);
        if isfinite(aliasFs)
            toneFreq = aliasToNyquist(tf(i), aliasFs);
        end
        [harmPow, ~, iHarm, iLeft, iRight] = signal.internal.getToneFromPSD(Pxx_bump_dc, F, rbw, toneFreq);  %%%
        [psdHarmPow(i), psdHarmFreq(i)] = idx2psddb(Pxx, F, iHarm);
        % obtain local maximum value in neighborhood of bin
        if ~isnan(harmPow)
            if ~isempty(iLeft) && ~isempty(iRight)
                % remove the power of this tone
                Pxx(iLeft(1):iRight(1)) = 0;
                harmIdx(i, :) = [iLeft(1); iRight(1)];
            end
        end
    end % of for i = 2:numel(tf)
end % of if length(tf) > 1
%%%%%%%%%%%%%%%%

% get an estimate of the noise floor by computing the median
% noise power of the non-harmonic region
estimatedNoiseDensity = median(Pxx(Pxx>0));

% extrapolate estimated noise density into dc/signal/harmonic regions
Pxx(Pxx==0) = estimatedNoiseDensity;

% prevent estimate from obscuring low peaks
Pxx = min([Pxx origPxx],[],2);

% compute the noise distortion.
totalNoise = bandpower(Pxx, F, 'psd');

r = 10*log10(Pfund / totalNoise);
noisePow = 10*log10(totalNoise);

if ~strcmp(plotType,'none')
    coder.internal.assert(coder.target('MATLAB'),'signal:snr:PlottingNotSupported');
    plotSNR(origPxx, F, rbw, plotType, psdHarmFreq, psdHarmPow, dcIdx, harmIdx);
    title(getString(message('signal:snr:SNRResult',sprintf('%6.2f',r))));
end

end

function toneF = aliasToNyquist(f, fs)

toneF1 = mod(f,fs);
if toneF1 > fs/2
    toneF = fs-toneF1;
else
    toneF = toneF1;
end

end

function plotSNR(Pxx, F, rbw, plotType, psdHarmFreq, psdHarmPow, dcIdx, harmIdx)

% scale Pxx by rbw
Pxx = Pxx * rbw;
psdHarmPow = psdHarmPow + 10*log10(rbw);

% initialize distortion plot
[hAxes, F, fscale, colors] = initdistplot(plotType, F);

% --- plot legend entry items ---

% plot fundamental
xData = F(harmIdx(1,1):harmIdx(1,2));
yData = 10*log10(Pxx(harmIdx(1,1):harmIdx(1,2)));
line(xData, yData, 'Parent', hAxes, 'Color', colors(1,:), 'LineWidth', 2);

% plot noise line
xData = F;
yData = 10*log10(Pxx);
line(xData, yData, 'Parent', hAxes, 'Color', colors(2,:));

% plot dc
xData = F(dcIdx(1):dcIdx(2));
yData = 10*log10(Pxx(dcIdx(1):dcIdx(2)));
line(xData, yData, 'Parent', hAxes, 'Color', colors(3,:), 'LineWidth', 2);

% --- use a solid grid slightly offset to accommodate text labels ---
initdistgrid(hAxes);

% --- replot over the grid ---
% plot dc marker
% mid_dcIdx = int64(floor((dcIdx(1)+dcIdx(2))/2));
% xData = F(mid_dcIdx)*fscale;
% yData = 10*log10(Pxx(mid_dcIdx));
% text(double(xData(1)),double(yData(1)),'DC', ...
%     'VerticalAlignment','bottom', ...
%     'HorizontalAlignment','center', ...
%     'BackgroundColor', 'w', ...
%     'EdgeColor', 'k', ...
%     'Color', colors(3,:));


% plot fundamental marker
xData = psdHarmFreq(1)*fscale;
yData = psdHarmPow(1);
text(double(xData(1)),double(yData(1)),'F', ...
    'VerticalAlignment','bottom', ...
    'HorizontalAlignment','center', ...
    'BackgroundColor', 'w', ...
    'EdgeColor', 'k', ...
    'Color', colors(1,:));

% plot harmonic markers
xData = double(psdHarmFreq(2:end)*fscale);
yData = double(psdHarmPow(2:end));
for i=1:numel(xData)    
    text(xData(i),yData(i),num2str(i+1), ...
        'VerticalAlignment','bottom', ...
        'HorizontalAlignment','center', ...
        'BackgroundColor', 'w', ...
        'EdgeColor', 'k', ...
        'Color', colors(3,:));
end

% plot noise line
xData = F;
yData = 10*log10(Pxx);
line(xData, yData, 'Parent', hAxes, 'Color', colors(2,:));

% plot fundamental
xData = F(harmIdx(1,1):harmIdx(1,2));
yData = 10*log10(Pxx(harmIdx(1,1):harmIdx(1,2)));
line(xData, yData, 'Parent', hAxes, 'Color', colors(1,:));

% plot dc on top
xData = F(dcIdx(1):dcIdx(2));
yData = 10*log10(Pxx(dcIdx(1):dcIdx(2)));
line(xData, yData, 'Parent', hAxes,'Color', colors(3,:));

% plot harmonics
for i=2:size(harmIdx,1)
    if ~any(isnan(harmIdx(i,:)))
        xData = F(harmIdx(i,1):harmIdx(i,2));
        yData = 10*log10(Pxx(harmIdx(i,1):harmIdx(i,2)));
        line(xData, yData, 'Parent', hAxes,'Color', colors(3,:));
    end
end

legend(getString(message('signal:snr:Fundamental')), ...
    getString(message('signal:snr:Noise')), ...
    getString(message('signal:snr:DCHarmonics')), ...
    'Location','best');

end


function [opt, idx] = localInputParser(validOpts,defaultOpt,arglist)

% Ensure any specified input string is a compile-time constant for code
% generation
coder.unroll();
for t = 1:numel(arglist)
    if ischar(arglist{t})
        coder.internal.assert(coder.internal.isConst(arglist{t}),...
            'signal:snr:FlagAsConst');
    end
end

found = false;
idx1 = 0;

coder.unroll();
for i = 1:numel(arglist)
    arg = arglist{i};
    coder.unroll();
    for j = 1:numel(validOpts)
        if ischar(arg) && strncmpi(arg,validOpts{j},length(arg))
            if ~found
                found = true;
                opt1 = validOpts{j};
                idx1 = i;
                break;
            end
            coder.internal.assert(~found,...
                'signal:snr:ConflictingOptions',opt1,validOpts{j});
        end
    end
end

if ~found
    opt = defaultOpt;
    idx = 0;
else
    opt = opt1;
    idx = idx1;
end

end

% LocalWords:  Bc Fs Pxx Sxx RBW NOISEPOW undersampled omitaliases Tpulse
% LocalWords:  sineex rbw enbw SINAD THD TOI undersampling replot

function chknostropts(varargin)
    %CHKNOSTROPTS Checks that no strings are specified
    %   errors on the first encountered string
    %
    %   This file is for internal use only
    
    %   Copyright 2013-2019 The MathWorks, Inc.
    
    %#codegen
    
    for i = 1:numel(varargin)
        coder.internal.errorIf(ischar(varargin{i}),...
            'signal:chknostropts:UnknownOption',varargin{i});
    end
end


function plotType = distplottype(n, esttype)
    %DISTPLOTTYPE Helper function for plotting power and psd estimates of
    %   distortion functions
    
    %   Copyright 2013-2019 The MathWorks, Inc.
    
    %#codegen
    
    if n>0
        plotType = 'none';
    elseif strcmp(esttype,'psd')
        plotType = 'psd';
    else
        plotType = 'power';
    end
end


function [opt, arglist] = getmutexclopt(validopts,defaultopt,arglist)
    %GETMUTEXCLOPT - get any of the specified mutually exclusive options and
    % remove from the argument list.  Allows initial matches.
    %
    %   This function is for internal purposes only and may be removed in a
    %   future release.
    %
    %   validtypes  - a cell array of valid options
    %                 (e.g. {'power','ms','psd'})
    %
    %   defaulttype - the default option to use if no type is found
    %
    %   arglist    - the input argument list
    %
    %   Errors out if different estimation types are matched in the arglist.    
    %
    %   See also CHKUNUSEDOPT.
    
    %   Copyright 2015 The MathWorks, Inc.
    
    opt = defaultopt;
    found = false;
    
    iarg = 1;
    while iarg <= numel(arglist)
      arg = arglist{iarg};
      if ischar(arg) && isrow(arg)
        matches = find(strncmpi(arg,validopts,length(arg)));
        if ~isempty(matches)
          if ~found
            found = true;
            opt = validopts{matches(1)};
            arglist(iarg) = [];
          else
            error(message('signal:getmutexclopt:ConflictingOptions', ...
                          opt,validopts{matches(1)}));
          end
        else
          iarg = iarg + 1;
        end
      else
       iarg = iarg + 1;
      end
    end
end

function [psd, freq] = idx2psddb(Pxx, F, idx)
    %IDX2PSDDB report back 10*log10(Pxx) and frequency at specified index
    %   return NaN when index is empty
    %
    %   This function is for internal use only. It may be removed in the future.
    
    %   Copyright 2013-2019 The MathWorks, Inc.
    
    %#codegen
    
    if ~isempty(idx)
        psd = 10*log10(Pxx(idx(1)));
        freq = F(idx(1));
    else
        psd = NaN('like',Pxx);
        freq = NaN('like',Pxx);
    end
end


function initdistgrid(hAxes)
    %INITDISTGRID initialize distortion plot grid
    %   Initialize the distortion grid so that it uses a light gray grid.
    %   Additionally make room for markers.
    %   
    %   This function is for internal use only. It may be removed in the future.
       
    %   Copyright 2013-2014 The MathWorks, Inc.
    
    % bound plot by y limit; give 10 extra dB for fundamental marker.
    yTick = get(hAxes,'YTick');
    if numel(yTick)>1
      yNew = yTick(end) + 10;
      yLim = get(hAxes,'YLim');
      set(hAxes,'YLim',[yLim(1) yNew]);
    end
    
    % Ensure axes limits are properly cached for zoom/unzoom
    resetplotview(hAxes,'SaveCurrentView');  
    
    % turn on box and grid
    set(hAxes,'Box','on','XGrid','on','YGrid','on');
end


function [hAxes, F, fscale, colorOrder] = initdistplot(plotType, F)
    %INITDISTPLOT initialize a plot for distortion annotation
    %
    %   This function is for internal use only. It may be removed in the future.
       
    %   Copyright 2013 The MathWorks, Inc.
    %%%%%%%%%%%%%
    figure;
    %%%%%%%%%%%%%
    hAxes = newplot;
    
    % set up x axis
    [F, fscale, xunits] = engunits(F);
    xlbl = getfreqlbl([xunits 'Hz']);
    xlabel(hAxes, xlbl);
    
    % bound plot by x limit
    set(hAxes, 'XLim', [min(F) max(F)]);
    set(hAxes, 'XScale', 'log'); %%%
    
    % set up y axis
    if strcmp(plotType,'power')
      ylbl = getString(message('signal:dspdata:dspdata:PowerdB'));
    else %'psd'
      ylbl = getString(message('signal:dspdata:dspdata:PowerfrequencydBHz'));
    end
    ylabel(hAxes, ylbl);
    
    % use the default color order (but set third color to blue instead of red)
    colorOrder = get(0,'DefaultAxesColorOrder');
    colorOrder(3,:) = [0 0 1];  % harmonic

    % fundamental: red
    colorOrder(1,:) = [1 0 0];

    % noise: 
    colorOrder(2,:) = [0 0 0];
end
