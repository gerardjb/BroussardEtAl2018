function [Decay_cf] = DecayFit_NonZero(Data,time_vect,plotOn,maskOut)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Joey Broussard
%UC Davis
%04/26/2018
%
%Function for fitting the decay kinetics of time series data. Data is fit
%by a double exponential and the plateau value is obligate as 0. The user
%can define a region to exclude from the fit if desired. Data is not
%smoothed.
%
%   Fit initialization is calculated as:
%       Initial value - average of the first three points
%       Plateau - average of the last three points
%       fast_tau - point at which data drops by 35%
%       slow_tau - point at which data drops by 85%. If tau points lie beyond
%           edge of data, data end point selected
%       slow/fast fraction - 50%
%
%   Inputs:
%       Data - time series data
%       time_vect - time vector. Default is indices of Data
%       plotOn - produce quality control plots. Default is off
%       maskOut - two element vector with first and last time points to
%           exclude. Default is no exclusion
%
%   Outputs:
%       tau_half - calculated time at which time series relaxes to 1/2 its
%           initial value
%       Decay_cf - decay curve fit object
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Setting up nargin conditions
if nargin==1
    time_vect = 1:length(Data);%indices if no time vector provided
    plotOn = 0;
    maskOut = [];
elseif nargin==2
    plotOn = 0;
    maskOut = [];
elseif nargin==3
    maskOut = [];
end

%Vertical arrange variables
time_vect = time_vect(:);
Data = Data(:);

%Establishing reasonable Start points for parameters
%End plateau is average of last three points
st_plateau = mean(Data((end-3):end));
%Establish reasonable guesses for taus if 50% frac is correct. This
%corresponds to intializing at a 35% and 85% drop, respectively.
st_init = mean(Data(1:3));
tau_fast_st_T = -0.35*(st_init) + st_init;
tau_slow_st_T = -0.85*(st_init) + st_init;
tau_fast_st = time_vect(find(Data<tau_fast_st_T,1));
%If trace decays too slow or too fast, set guess to last or second time
%point, respectively
if isempty(tau_fast_st)
    tau_fast_st = time_vect(end);
elseif tau_fast_st==0
    tau_fast_st = time_vect(2);
end
tau_slow_st = time_vect(find(Data<tau_slow_st_T,1));
if isempty(tau_slow_st)
    tau_slow_st = time_vect(end);
elseif tau_slow_st==0
    tau_slow_st = time_vect(2);
end

%Setting limits and initial conditions for parameters
fo = fitoptions('method','NonlinearLeastSquares',...
    'Lower',[-10000,-10000,0,0,0],...
    'Upper',[10000,10000,1,1000,2000]);
ok_ = isfinite(time_vect) & isfinite(Data);
%Excluded regions
if ~isempty(maskOut)
    %setting excluded points
    exclude = false(length(time_vect),1);
    exclude(time_vect>maskOut(1) & time_vect<maskOut(2)) = 1;
    %initializing excluded points
    set(fo,'Exclude',exclude(ok_));
end
%Initialize parameter values
st = [st_init,st_plateau,0.5,tau_fast_st,tau_slow_st];
set(fo,'Startpoint',st);

%model to fit to
ft = fittype(...
    'init - (init - plateau)*(1-FastFrac*exp(-t/tau_fast)-(1-FastFrac)*exp(-t/tau_slow))',...
    'dependent',{'y'},'independent',{'t'},...
    'coefficients',{'init', 'plateau', 'FastFrac', 'tau_fast', 'tau_slow'});

%Fit the model using selected data
Decay_cf = fit(time_vect,Data,ft,fo);


%Make qc plot if user desires
if plotOn
    %get plottable version of data
    decay_hat = feval(Decay_cf,time_vect);
    
    %plot of data with fit
    plot(time_vect,Data,'-r.',...
        time_vect,decay_hat,'k');
    title('Decay curve fit to time series data')
    ylabel('time series data')
    xlabel('time (A.U.)')
    legend('data series','double exponential fit')
end

end

