%% HHM
% Hybrid Harmonic Model (HHM)
%
% Purpose
%   This script builds long time series of total water level (TWL) for a chosen
%   location by combining:
%     1) Astronomical tides reconstructed from harmonic constituents (t tide)
%     2) Mean sea level change (SLR scenario curve)
%     3) Storm surge and non tidal residual (NTR) from Monte Carlo simulations
%     4) Optional wind wave setup and runup (if computewindwaves equals 1)
%
% Main outputs (variables created in the MATLAB workspace)
%   WL         Astronomical water level plus SLR and datum shift, size Nt by 3
%              per rate. Three columns represent low, median, high SLR paths
%   Twl_r20    TWL cell array, 1 by 3, each cell is Nt by Nsims
%   Twl_r37    TWL cell array, 1 by 3, each cell is Nt by Nsims
%   Twl_r52    TWL cell array, 1 by 3, each cell is Nt by Nsims
%   Twl_r100   TWL cell array, 1 by 3, each cell is Nt by Nsims
%   Optional: RS_timeseries is created when computewindwaves equals 1
%
% Input data used by this script
%   A) Delft3D FM history NetCDF files, variable name waterlevel
%      File pattern is built in the loop under the section titled
%      MODEL WATER LEVELS RATES times SLR paths
%      Each file is read using ncread and then subset to a single mesh index
%   B) Storm surge Monte Carlo MAT files
%      all_montecarlo_SS.mat must contain ss_data_matrix
%      montecarlo_time.mat must contain time_ss
%   C) Wind wave setup and runup lookup tables when computewindwaves equals 1
%      RSmat_1_2m_newbathy_v2.mat and RSmat_2_3point3m_newbathy_v2.mat
%
% External dependencies
%   t tide toolbox (t_tide and t_predic) must be on the MATLAB path
%   Inpaint_nans is used by the wind wave functions path block
%   get_runupsetup_timeseries_v6 must be on the MATLAB path when wind waves run
%
% How to run
%   1) Set the user settings in the next section
%   2) Ensure the input folders and MAT files exist
%   3) Run the script
%
% Notes on phase interpolation
%   Constituent phase is unwrapped in radians, interpolated with PCHIP, and
%   wrapped back to 0 to 360 degrees. This preserves phase continuity.

clc; clear all; close all;

%% User settings
% Edit only this block for a clean GitHub release
computewindwaves = 1; % set to 0 to skip wind wave setup and runup

% Location of the t tide toolbox folder
% Example: 'D:\\RESEARCH\\HARMONIC_MODEL\\t_tide'
t_tide_folder = 'D:\\RESEARCH\\HARMONIC_MODEL\\t_tide';

% Working directory where the relative paths in this script are valid
% Example: 'D:\\RESEARCH\\HARMONIC_MODEL\\Model_Files\\Model_Files'
project_folder = 'D:\\RESEARCH\\HARMONIC_MODEL\\Model_Files\\Model_Files';

% Optional output save
save_outputs = 1;               % set to 0 to disable saving to MAT
output_folder = fullfile(project_folder, 'outputs');

% Datum shift applied in the tidal reconstruction (MSL to NAVD88)
msl_to_navd88_m = 1.18;

% Basic input checks with clear error messages
if ~exist(t_tide_folder,'dir')
  error('t_tide_folder not found. Update t_tide_folder in the User settings block.')
end
if ~exist(project_folder,'dir')
  error('project_folder not found. Update project_folder in the User settings block.')
end

%% Parallel setup that never overrequests workers and falls back cleanly
% Tries Processes cluster then local
% Uses up to 8 workers but respects profile limit and CPU cores
% Falls back to serial if Parallel Computing Toolbox is unavailable
try
  try
    c = parcluster('Processes');
  catch
    c = parcluster('local');
  end
  maxAllowed = c.NumWorkers;
  nw = min([4, feature('numcores'), maxAllowed]);
  p = gcp('nocreate');
  if isempty(p) || p.NumWorkers ~= nw
    parpool(c, nw);
  end
catch ME
  warning(['Parallel pool not started. Running serially. Details: ' ME.message]);
end

%% Define time, load datasets and add t tide path
% Add required toolboxes and move to the project working folder
addpath(t_tide_folder)
cd(project_folder)

t = datetime(2021,12,1):minutes(10):datetime(2022,1,1);
% NOTE NEW MODEL
% US101 equals 95
% AIRPORT equals 31
% DOWNTOWN equals 74
% AGRIC AREA equals 92 or 220 depending on your mesh reference

%% ================== DYNAMIC LOCATION CONTROL ==================
% Single switch for the analysis location. Everything else updates from this choice
% location_waves equals 1 means downtown, 2 means 101, 3 means airport, 4 means agricultural area
location_waves = 3; % 1 downtown, 2 US-101, 3 airport, 4 agricultural area

% Map location choice to idx, titles, legends, and threshold
switch location_waves
  case 1 % downtown
    idx = 74;
    loc_title = 'Downtown';
    legend_pd = 'Present Day Levee Downtown';
    legend_fwr = 'Full Wetland Restoration Downtown';
    exceed_threshold_m = 3.25; % DT equals 3.25 m NAVD88
  case 2 % US 101
    idx = 95;
    loc_title = 'US 101 Highway';
    legend_pd = 'Present Day Levee US 101 Highway';
    legend_fwr = 'Full Wetland Restoration US 101 Highway';
    exceed_threshold_m = 3.32; % Highway equals 3.32 m NAVD88
  case 3 % Airport
    idx = 31;
    loc_title = 'Airport';
    legend_pd = 'Present Day Levee Airport';
    legend_fwr = 'Full Wetland Restoration';
    exceed_threshold_m = 4; % set your site specific value if different
  case 4 % Agricultural area
    idx = 220; % confirmed agric index is 220
    loc_title = 'Agric';
    legend_pd = 'Present Day Levee Agric';
    legend_fwr = 'Full Wetland Restoration Agric';
    exceed_threshold_m = 3.25;
  otherwise
    error('location_waves must be 1, 2, 3, or 4');
end

%% ================== MODEL WATER LEVELS RATES times SLR paths ==================
rates  = {'r20','r37','r52','r100'};
slrtag = {'l','m','h'};  % low, median, high SLR
years5 = [2020 2040 2060 2080 2100];
slrFolder = struct('l','L','m','M','h','H');

for ir = 1:numel(rates)
  rate = rates{ir};
  for is = 1:numel(slrtag)
    s = slrtag{is};
    dataCols = cell(1,5);
    for iy = 1:5
      yr = years5(iy);
      if yr == 2020
        fdir = sprintf('./GMSLR_1m/s%04d_%s/output/FlowFM_0000_his.nc', yr, rate);
      else
        fdir = sprintf('./GMSLR_1m/s%04d%s_%s/output/FlowFM_0000_his.nc', yr, slrFolder.(s), rate);
      end
      dum = ncread(fdir,'waterlevel');
      dataCols{iy} = dum(idx,1153:end)'; % December start index
    end
    R.(rate).(s) = [dataCols{:}];
  end
end

%% Calculate SLR fits
t = 1:81;
SLR.scenario = 1;
SLR.fit = fit([0 20 40 60 80]',[0 0.08 0.18 0.36 0.62 ]','poly4');
SLR.msl_1 = SLR.fit(t);
SLR.scenario = 2;
SLR.fit = fit([0 20 40 60 80]',[0 0.11 0.24 0.45 0.81 ]','poly4');
SLR.msl_2 = SLR.fit(t);
SLR.scenario = 3;
SLR.fit = fit([0 20 40 60 80]',[0 0.16 0.34 0.58 0.93 ]','poly4');
SLR.msl_3 = SLR.fit(t);
SLR.msl = [SLR.msl_1, SLR.msl_2, SLR.msl_3];

%% ================== Calculate changes to tidal constituents over time ==================
x_years = years5; Ny = numel(x_years);
for ir = 1:numel(rates)
  rate = rates{ir};
  tidestruc_c = cell(1,Ny);
  parfor k = 1:Ny
    tidestruc_c{k} = t_tide(R.(rate).l(:,k), 'interval',1/6,'error','wboot','output','none','latitude',43.35,...
                            'start time',datenum(x_years(k),12,1));
  end
  T.(rate).l = [tidestruc_c{:}];
  tidestruc_c = cell(1,Ny);
  parfor k = 1:Ny
    tidestruc_c{k} = t_tide(R.(rate).m(:,k), 'interval',1/6,'error','wboot','output','none','latitude',43.35,...
                            'start time',datenum(x_years(k),12,1));
  end
  T.(rate).m = [tidestruc_c{:}];
  tidestruc_c = cell(1,Ny);
  parfor k = 1:Ny
    tidestruc_c{k} = t_tide(R.(rate).h(:,k), 'interval',1/6,'error','wboot','output','none','latitude',43.35,...
                            'start time',datenum(x_years(k),12,1));
  end
  T.(rate).h = [tidestruc_c{:}];
end

cd(project_folder)
%% ================== Interpolate between years (PHASE-CONTINUOUS PCHIP, NO DIPS) ==================
% UPDATED FIX:
%   Replaces the complex-plane interpolation with separate amplitude and phase interpolation
%   using unwrap() to preserve phase continuity and avoid unphysical amplitude dips caused
%   by vector cancellation.
%
%   Procedure:
%     1) Unwrap phase in radians (shortest rotation)
%     2) Interpolate amplitude and unwrapped phase separately using PCHIP
%     3) Reconvert phase to degrees and wrap to [0..360]
%   Uncertainties (ampu, phaseu) still interpolated directly with PCHIP.

x  = [2020, 2040, 2060, 2080, 2100];
xq = 2020:2100;

for ir = 1:numel(rates)
  rate = rates{ir};

  % Names/freq/type for t_predic
  ftidestruc.(rate).name = T.(rate).h(1).name;
  ftidestruc.(rate).freq = T.(rate).h(1).freq;
  ftidestruc.(rate).type = 'nodal';

  % Gather the 5-year tidecon slabs: [constituent, (amp ampu phase phaseu), year, SLR]
  tidestruc5 = zeros(29,4,5,3);
  for k = 1:5, tidestruc5(:,:,k,1) = T.(rate).l(k).tidecon; end
  for k = 1:5, tidestruc5(:,:,k,2) = T.(rate).m(k).tidecon; end
  for k = 1:5, tidestruc5(:,:,k,3) = T.(rate).h(k).tidecon; end

  % Output array for 2020..2100 (81 years) per SLR branch
  ftidestruc.(rate).tidecon = zeros(29,4,81,3);

  for s = 1:3
    for i = 1:29
      % 5-year amplitude and phase (degrees)
      amp5 = squeeze(tidestruc5(i,1,:,s));      % amplitude >= 0
      ph5  = squeeze(tidestruc5(i,3,:,s));      % phase in degrees [0..360]

      % --- Fix: unwrap phase and interpolate separately ---
      ph5u = unwrap(deg2rad(ph5));              % unwrap once (shortest direction)
      ampq = pchip(x, amp5, xq);                % amplitude interpolation (smooth, no overshoot)
      phq  = pchip(x, ph5u, xq);                % phase interpolation in radians

      % Convert back to degrees and wrap to [0,360)
      ftidestruc.(rate).tidecon(i,1,:,s) = ampq;
      ftidestruc.(rate).tidecon(i,3,:,s) = mod(rad2deg(phq), 360);

      % Interpolate uncertainties directly with PCHIP
      ftidestruc.(rate).tidecon(i,2,:,s) = pchip(x, squeeze(tidestruc5(i,2,:,s)), xq); % ampu
      ftidestruc.(rate).tidecon(i,4,:,s) = pchip(x, squeeze(tidestruc5(i,4,:,s)), xq); % phaseu
    end
  end
end

%% Project future tidal water levels using interpolated constituents
% t100 is datetime 2020 01 01 through 2101 01 01 with 10 minute step
t100 = [datetime(2020,1,1):minutes(10):datetime(2101,1,1)]';
for ir = 1:numel(rates)
  rate = rates{ir};
  WL.(rate) = zeros([length(t100),3]);
  c = 1;
  for j = 1:3 % SLR paths
    for i = 1:length(xq)-1
      dum = t_predic( datenum(xq(i),1,1):0.0069444444444444:datenum(xq(i+1),1,1), ...
                      ftidestruc.(rate).name, ftidestruc.(rate).freq, ftidestruc.(rate).tidecon(:,:,i,j), ...
                      'latitude', 43.35 );
      tlen = length(dum);
      WL.(rate)(c:c+tlen-1,j) = dum + SLR.msl(i,j) + msl_to_navd88_m; % WL plus SLR scenario plus datum shift
      c = c + tlen;
    end
    c = 1;
  end
end

%% ================== LOAD STORM SURGE DATA ==================
if ~exist('all_montecarlo_SS.mat','file')
  error('Missing input MAT file all_montecarlo_SS.mat. Place it in project_folder or update the load path.')
end
if ~exist('montecarlo_time.mat','file')
  error('Missing input MAT file montecarlo_time.mat. Place it in project_folder or update the load path.')
end

load('all_montecarlo_SS.mat'); % 100 Monte Carlo simulations 2020 through 2120
NTR_2020_2100 = ss_data_matrix(1:701281, :); % extract data matching 2020 through 2100 hourly
load('montecarlo_time.mat'); % time matrix for 100 Monte Carlo simulations

time_ss = time_ss(1:701281, :);
time_reference = '1970-01-01 00:00:00'; % time in datetime
NTR_t100 = datetime(time_reference, 'InputFormat', 'yyyy-MM-dd HH:mm:ss') + days(time_ss);

% IMPORTANT: restore your original behavior here (linear interpolation)
NTR_interpolated = interp1(NTR_t100, NTR_2020_2100, t100, 'linear');   % <- back to linear


%% ================== wind calculation here Sam ==================
% load water level and depth as raster if needed
% grid_res equals 5 m
% [Zq, Zq_bl] equals loadtopobathy grid_res
addpath(fullfile(project_folder,'Inpaint_nans','Inpaint_nans'))
if computewindwaves==1
  if ~exist('RSmat_1_2m_newbathy_v2.mat','file')
    error('Missing input MAT file RSmat_1_2m_newbathy_v2.mat. Place it in project_folder or update the load path.')
  end
  if ~exist('RSmat_2_3point3m_newbathy_v2.mat','file')
    error('Missing input MAT file RSmat_2_3point3m_newbathy_v2.mat. Place it in project_folder or update the load path.')
  end

  % DRIVER_windwaves_coos_v5 prepares matrices used below
  load RSmat_1_2m_newbathy_v2.mat
  Rmat1_2 = Rmat(1:11,:,:,:)./1.29; % convert to R2 percent instead of Rmax
  Smat1_2 = Smat(1:11,:,:,:);
  load RSmat_2_3point3m_newbathy_v2.mat
  Rmat2_3p3 = Rmat(1:13,:,:,:)./1.29;
  Smat2_3p3 = Smat(1:13,:,:,:);
  Rmatfull  = cat(1, Rmat1_2,  Rmat2_3p3);
  Smatfull  = cat(1, Smat1_2,  Smat2_3p3);
  if     idx==74  % downtown
      Rmatproc = Rmatfull(:,:,:,1);  Smatproc = Smatfull(:,:,:,1);
  elseif idx==95  % US101
      Rmatproc = Rmatfull(:,:,:,2);  Smatproc = Smatfull(:,:,:,2);
  elseif idx==31  % airport
      Rmatproc = Rmatfull(:,:,:,3);  Smatproc = Smatfull(:,:,:,3);
  elseif idx==220 % agric
      Rmatproc = Rmatfull(:,:,:,4);  Smatproc = Smatfull(:,:,:,4);
  end
  clear Rmatfull Smatfull

  % Parallel Monte Carlo loop with Constant blocks to minimize broadcasting overhead
  Nsims = size(NTR_interpolated,2);
  Nt    = size(WL.(rates{1}),1);

  % dimensions tstep by locations DT 101 APT AGR by scenario L M H by simulation
  for ir = 1:numel(rates)
    rate = rates{ir};
    RS_timeseries.(rate) = zeros(Nt, 3, 1, Nsims, 'like', WL.(rate));
  end

  cNTR = parallel.pool.Constant(NTR_interpolated);
  cR   = parallel.pool.Constant(Rmatproc);
  cS   = parallel.pool.Constant(Smatproc);
  cT   = parallel.pool.Constant(t100);

  for ir = 1:numel(rates)
    rate = rates{ir};
    cWL  = parallel.pool.Constant(WL.(rate));
    RS_tmp = zeros(Nt,3,1,Nsims,'like',WL.(rate));
    parfor numsim = 1:Nsims
        randsim = numsim; % same as number of Monte Carlos
        WLsim   = cWL.Value + cNTR.Value(:,numsim); % WL for waves is HHM plus NTR storm surges all SLR at once
        RS_sim  = get_runupsetup_timeseries_v6(WLsim, cT.Value, cR.Value, cS.Value, randsim); % dimensions tstep by locations by scenario L M H
        RS_tmp(:,:,:,numsim) = RS_sim;
    end
    RS_tmp(isnan(RS_tmp)==1) = 0; % set infrequent NaN to zero
    RS_timeseries.(rate) = RS_tmp;
  end
end
% end Sam

%% ================== Calc astronomic water level plus NTR plus optional runup for each SLR scenario ==================
% Creates scenario specific TWL outputs that carry the scenario number
% Twl_r20, Twl_r37, Twl_r52, Twl_r100
% Each is a 1 by 3 cell array Low SLR, Med SLR, High SLR and each cell is Nt by Nsim
for ir = 1:numel(rates)
  rate = rates{ir};
  [~, num_scenarios] = size(WL.(rate));  %#ok<ASGLU>
  results_cell = cell(1,num_scenarios);

  for i = 1:num_scenarios
      base_wl = WL.(rate)(:, i); % Nt by 1
      if computewindwaves==1
          RS_loc = squeeze(RS_timeseries.(rate)(:, :, 1, :));  % Nt by 3 by Nsims
          results_cell{i} = base_wl + NTR_interpolated + squeeze(RS_loc(:, i, :));
      else
          % RUNUP OFF so TWL equals WL plus NTR only
          results_cell{i} = base_wl + NTR_interpolated;
      end
  end

  switch rate
    case 'r20',  Twl_r20  = results_cell; %#ok<NASGU>
    case 'r37',  Twl_r37  = results_cell; %#ok<NASGU>
    case 'r52',  Twl_r52  = results_cell; %#ok<NASGU>
    case 'r100', Twl_r100 = results_cell; %#ok<NASGU>
  end
end
        
%% ================== Calculate hours of exceedances ==================
% IND equals 1 low, 2 med, 3 high
% Define the exceedance threshold from the location block
threshold = exceed_threshold_m;  % Highway equals 3.32 m and DT equals 3.25 m NAVD88
dt_hours  = 10/60;

% Helper to compute hours across simulations for a 1 by 3 TWL cell
compute_hours = @(C) vertcat( ...
  sum(C{1} > threshold, 1) * dt_hours, ...
  sum(C{2} > threshold, 1) * dt_hours, ...
  sum(C{3} > threshold, 1) * dt_hours);

ex.dt_hr_r20  = compute_hours(Twl_r20);
ex.dt_hr_r37  = compute_hours(Twl_r37);
ex.dt_hr_r52  = compute_hours(Twl_r52);
ex.dt_hr_r100 = compute_hours(Twl_r100);

%% ================== Plotting exceedances across the TWL Monte Carlo simulations ==================
% Bar chart of mean plus std hours of exceedance for three SLR paths and four rate bars
% Color mapping
% r52 blue Top 10 sites restored
% r100 green Full wetland restoration
% r37 purple Two largest sites restored
% r20 orange Present day
col.r20  = [1.000 0.549 0.000];  % orange
col.r37  = [0.494 0.184 0.556];  % purple
col.r52  = [0.150 0.370 0.730];  % blue
col.r100 = [0.000 0.550 0.200];  % green

% means and stds with omitnan
mean_dt = [ ...
  mean(ex.dt_hr_r20 ,2,'omitnan'), ...
  mean(ex.dt_hr_r37 ,2,'omitnan'), ...
  mean(ex.dt_hr_r52 ,2,'omitnan'), ...
  mean(ex.dt_hr_r100,2,'omitnan') ];

std_dt  = [ ...
  std(ex.dt_hr_r20 ,0,2,'omitnan'), ...
  std(ex.dt_hr_r37 ,0,2,'omitnan'), ...
  std(ex.dt_hr_r52 ,0,2,'omitnan'), ...
  std(ex.dt_hr_r100,0,2,'omitnan') ];

SLR_labels = {'SLR = 62 cm','SLR = 81 cm','SLR = 93 cm'};

figure
b = bar(mean_dt,'grouped'); hold on
b(1).FaceColor = col.r20;   b(1).DisplayName = 'R20 Present day';
b(2).FaceColor = col.r37;   b(2).DisplayName = 'R37 Two largest restoration';
b(3).FaceColor = col.r52;   b(3).DisplayName = 'R52 Top 10 restoration';
b(4).FaceColor = col.r100;  b(4).DisplayName = 'R100 Full restoration';

% real error bars: include ONLY the first one in the legend
hErr = gobjects(1,numel(b));
for ib = 1:numel(b)
  hErr(ib) = errorbar(b(ib).XEndPoints, mean_dt(:,ib), std_dt(:,ib), ...
                      'k', 'linestyle','none');
  if ib == 1
    set(hErr(ib),'DisplayName','Standard Deviation');  % legend entry
  else
    set(hErr(ib),'HandleVisibility','off');            % keep legend clean
  end
end

ylabel('Hours of Flooding','FontWeight','bold');
set(gca,'XTickLabel',SLR_labels);

legend([b hErr(1)], {'R20 Present day','R37 Two largest restoration', ...
                     'R52 Top 10 restoration','R100 Full restoration', ...
                     'Standard Deviation'}, 'Location','northwest');

title(['Hours of Levee Exceedance ' (loc_title)]);
ax = gca; ax.FontSize = 12; grid on; hold off

%% Decadal Exceedance Plot
%% Decadal exceedance one figure all rates with SLR scenarios overlaid

% use these inputs that already exist in your workspace
% Twl_r20 Twl_r37 Twl_r52 Twl_r100  each is 1 by 3 cell array (low med high)
% t100                                datetime vector
% exceed_threshold_m                  levee threshold (m)
% col.r20 col.r37 col.r52 col.r100    colour map from earlier

threshold = exceed_threshold_m;
dt_hours  = 10/60;          % 10 minute model time step in hours

% custom decades 2040 through 2100
decade_edges  = [2040 2050 2060 2070 2080 2090 2100];
D             = numel(decade_edges) - 1;
decade_labels = string(decade_edges(2:end));   % 2050 2060 2070 2080 2090 2100

rates     = {'r20','r37','r52','r100'};
Twlcells  = {Twl_r20, Twl_r37, Twl_r52, Twl_r100};
col_rate  = {col.r20, col.r37, col.r52, col.r100};

% mean_decade(ir, scen, d)
% error values stored as distance from mean to min and mean to max
mean_decade     = zeros(4,3,D);
err_low_decade  = zeros(4,3,D);
err_high_decade = zeros(4,3,D);

for ir = 1:4                         % restoration rate index
  Acell = Twlcells{ir};              % cell with three SLR paths

  for s_ind = 1:3                    % 1 low 2 med 3 high SLR
    A = Acell{s_ind};                % Nt by Nsims matrix

    for d = 1:D
      % decade window
      t0   = datetime(decade_edges(d),   1, 1);
      t1   = datetime(decade_edges(d+1), 1, 1) + caldays(-1);
      idxD = (t100 >= t0 & t100 <= t1);

      % hours of exceedance for each Monte Carlo run
      hrs_vec = sum(A(idxD,:) > threshold, 1) * dt_hours;   % 1 by Nsims

      m  = mean(hrs_vec,'omitnan');
      mn = min(hrs_vec,[],'omitnan');
      mx = max(hrs_vec,[],'omitnan');

      mean_decade    (ir,s_ind,d) = m;
      err_low_decade (ir,s_ind,d) = m - mn;
      err_high_decade(ir,s_ind,d) = mx - m;
    end
  end
end

xbase    = 1:D;                 % one group per decade
w        = 0.20;                % bar width

off = [-1.62 -0.54 0.54 1.62] * w; % offsets for R20 R37 R52 R100 within each decade

alphas   = [1.0 0.50 0.3];      % low med high SLR transparency

RateLabels = {'R20 Present day', ...
              'R37 Two largest restoration', ...
              'R52 Top 10 restoration', ...
              'R100 Full restoration'};

figure;
hold on;

% one handle per restoration rate for the legend
hRate = gobjects(4,1);

for d = 1:D                      % decade loop
  for ir = 1:4                   % restoration rate loop
    x0 = xbase(d) + off(ir);     % x position for this rate in this decade

    for s_ind = 1:3              % SLR branch loop
      mval  = mean_decade    (ir,s_ind,d);
      elow  = err_low_decade (ir,s_ind,d);
      ehigh = err_high_decade(ir,s_ind,d);

      h = bar(x0, mval, w, ...
        'FaceColor',  col_rate{ir}, ...
        'FaceAlpha',  alphas(s_ind), ...
        'EdgeColor',  'none');

      % asymmetric error bars min and max around the mean
      errorbar(x0, mval, elow, ehigh, 'k', 'LineStyle','none');

      % keep only low SLR bar from first decade for legend
      if d == 1 && s_ind == 1
        hRate(ir) = h;
      end

      % hide medium and high SLR bars in legend
      if s_ind ~= 1
        h.HandleVisibility = 'off';
      end
    end
  end
end

% axes formatting
xticks(xbase);
xticklabels(decade_labels);
xlabel('Decade End Year','FontWeight','bold');
ylabel('Hours of Flooding','FontWeight','bold');
title(['Levee Exceedance Per Decade ' (loc_title)],'FontWeight','bold');
grid on;
set(gca,'FontSize',12);

% force y axis to start at zero but leave top automatic
ax = gca;

yl = ax.YLim;
ax.YLim = [0 yl(2)];

% legend only four entries one per restoration rate
legend(hRate, RateLabels, 'Location','northwest');

hold off;

%% ================== Annual exceedance ==================

%% %% Annual hours of exceedance with SLR bands 2025-2100 (median across Monte Carlo)

% ---------- user controls for plot window and ticks ----------
year_start_plot = 2025;   % first year to show on x axis (you can set 2020, 2030, etc)
year_end_plot   = 2099;   % last year to show on x axis (for full record use 2100)
year_tick_step  = 5;      % x tick spacing in years
% -------------------------------------------------------------
y_start= 2025;
y_end =2100;
dt_hours = 10/60;
all_years   = year(t100(1)):year(t100(end)) - 1;  % complete years 2020 2021 ... 2100
num_years   = numel(all_years);

% containers: (slr_index, year)
ann_r20  = zeros(3, num_years);
ann_r37  = zeros(3, num_years);
ann_r52  = zeros(3, num_years);
ann_r100 = zeros(3, num_years);

for s_ind = 1:3                 % 1 low 2 med 3 high SLR

  A20  = Twl_r20{s_ind};        % Nt by Nsims
  A37  = Twl_r37{s_ind};
  A52  = Twl_r52{s_ind};
  A100 = Twl_r100{s_ind};

  for iy = 1:num_years
    yr_mask = (year(t100) == all_years(iy));

    hrs20_vec  = sum(A20 (yr_mask,:) > threshold, 1) * dt_hours;
    hrs37_vec  = sum(A37 (yr_mask,:) > threshold, 1) * dt_hours;
    hrs52_vec  = sum(A52 (yr_mask,:) > threshold, 1) * dt_hours;
    hrs100_vec = sum(A100(yr_mask,:) > threshold, 1) * dt_hours;

    % median Monte Carlo value for this SLR branch
    % if you want the mean instead, replace median with mean
    ann_r20 (s_ind,iy)  = median(hrs20_vec,  'omitnan');
    ann_r37 (s_ind,iy)  = median(hrs37_vec,  'omitnan');
    ann_r52 (s_ind,iy)  = median(hrs52_vec,  'omitnan');
    ann_r100(s_ind,iy)  = median(hrs100_vec, 'omitnan');
  end
end

% restrict to user specified plotting window
idx_plot    = (all_years >= year_start_plot) & (all_years <= year_end_plot);
years_plot  = all_years(idx_plot);

R20_low   = ann_r20(1,idx_plot);
R20_med   = ann_r20(2,idx_plot);
R20_high  = ann_r20(3,idx_plot);

R37_low   = ann_r37(1,idx_plot);
R37_med   = ann_r37(2,idx_plot);
R37_high  = ann_r37(3,idx_plot);

R52_low   = ann_r52(1,idx_plot);
R52_med   = ann_r52(2,idx_plot);
R52_high  = ann_r52(3,idx_plot);

R100_low  = ann_r100(1,idx_plot);
R100_med  = ann_r100(2,idx_plot);
R100_high = ann_r100(3,idx_plot);

%% plotting: low high band plus medium SLR line for each rate

figure; hold on;

x_band = [years_plot, fliplr(years_plot)];

% R20 band and medium line
y_lo = min(R20_low, R20_high);
y_hi = max(R20_low, R20_high);
hR20band = fill(x_band, [y_lo, fliplr(y_hi)], col.r20, ...
                'FaceAlpha',0.20, 'EdgeColor','none');
hR20med  = plot(years_plot, R20_med, 'Color',col.r20, 'LineWidth',2);

% R37 band and medium line
y_lo = min(R37_low, R37_high);
y_hi = max(R37_low, R37_high);
hR37band = fill(x_band, [y_lo, fliplr(y_hi)], col.r37, ...
                'FaceAlpha',0.20, 'EdgeColor','none');
hR37med  = plot(years_plot, R37_med, 'Color',col.r37, 'LineWidth',2);

% R52 band and medium line
y_lo = min(R52_low, R52_high);
y_hi = max(R52_low, R52_high);
hR52band = fill(x_band, [y_lo, fliplr(y_hi)], col.r52, ...
                'FaceAlpha',0.20, 'EdgeColor','none');
hR52med  = plot(years_plot, R52_med, 'Color',col.r52, 'LineWidth',2);

% R100 band and medium line
y_lo = min(R100_low, R100_high);
y_hi = max(R100_low, R100_high);
hR100band = fill(x_band, [y_lo, fliplr(y_hi)], col.r100, ...
                 'FaceAlpha',0.20, 'EdgeColor','none');
hR100med  = plot(years_plot, R100_med, 'Color',col.r100, 'LineWidth',2);

% make sure lines sit above shading
uistack([hR20med hR37med hR52med hR100med],'top');

grid on;
xlabel('Year','FontWeight','bold');
ylabel('Hours of Flooding','FontWeight','bold');
title(sprintf('Annual hours of exceedance %d to %d (%s)', ...
      y_start, y_end, loc_title), 'FontWeight','bold');

set(gca,'FontSize',12);

% y axis from zero upward, top auto
ax = gca;
yl = ax.YLim;
ax.YLim = [0 yl(2)];



% x axis limits and ticks tied to user settings
xlim([year_start_plot year_end_plot]);
xticks(year_start_plot:year_tick_step:year_end_plot);
xtickangle(45);
% legend formatted like your example figure
LG = legend([hR20band hR20med ...
        hR37band hR37med ...
        hR52band hR52med ...
        hR100band hR100med], ...
       {'R20 Low-High SLR band','R20 Med SLR', ...
        'R37 Low-High SLR band','R37 Med SLR', ...
        'R52 Low-High SLR band','R52 Med SLR', ...
        'R100 Low-High SLR band','R100 Med SLR'}, ...
       'Location','northwest', 'NumColumns',2);
set(LG,'FontSize',10);   % smaller text
hold off;

%% ANNUAL 2025-2050
% ---------- user controls for plot window and ticks ----------
year_start_plot = 2025;   % first year to show on x axis (you can set 2020, 2030, etc)
year_end_plot   = 2050;   % last year to show on x axis (for full record use 2100)
year_tick_step  = 5;      % x tick spacing in years
% -------------------------------------------------------------
y_startt= 2025;
y_endd =2050;
dt_hours = 10/60;
all_years   = year(t100(1)):year(t100(end)) - 1;  % complete years 2020 2021 ... 2100
num_years   = numel(all_years);

% containers: (slr_index, year)
ann_r20  = zeros(3, num_years);
ann_r37  = zeros(3, num_years);
ann_r52  = zeros(3, num_years);
ann_r100 = zeros(3, num_years);

for s_ind = 1:3                 % 1 low 2 med 3 high SLR

  A20  = Twl_r20{s_ind};        % Nt by Nsims
  A37  = Twl_r37{s_ind};
  A52  = Twl_r52{s_ind};
  A100 = Twl_r100{s_ind};

  for iy = 1:num_years
    yr_mask = (year(t100) == all_years(iy));

    hrs20_vec  = sum(A20 (yr_mask,:) > threshold, 1) * dt_hours;
    hrs37_vec  = sum(A37 (yr_mask,:) > threshold, 1) * dt_hours;
    hrs52_vec  = sum(A52 (yr_mask,:) > threshold, 1) * dt_hours;
    hrs100_vec = sum(A100(yr_mask,:) > threshold, 1) * dt_hours;

    % median Monte Carlo value for this SLR branch
    % if you want the mean instead, replace median with mean
    ann_r20 (s_ind,iy)  = median(hrs20_vec,  'omitnan');
    ann_r37 (s_ind,iy)  = median(hrs37_vec,  'omitnan');
    ann_r52 (s_ind,iy)  = median(hrs52_vec,  'omitnan');
    ann_r100(s_ind,iy)  = median(hrs100_vec, 'omitnan');
  end
end

% restrict to user specified plotting window
idx_plot    = (all_years >= year_start_plot) & (all_years <= year_end_plot);
years_plot  = all_years(idx_plot);

R20_low   = ann_r20(1,idx_plot);
R20_med   = ann_r20(2,idx_plot);
R20_high  = ann_r20(3,idx_plot);

R37_low   = ann_r37(1,idx_plot);
R37_med   = ann_r37(2,idx_plot);
R37_high  = ann_r37(3,idx_plot);

R52_low   = ann_r52(1,idx_plot);
R52_med   = ann_r52(2,idx_plot);
R52_high  = ann_r52(3,idx_plot);

R100_low  = ann_r100(1,idx_plot);
R100_med  = ann_r100(2,idx_plot);
R100_high = ann_r100(3,idx_plot);

%% plotting: low high band plus medium SLR line for each rate

figure; hold on;

x_band = [years_plot, fliplr(years_plot)];

% R20 band and medium line
y_lo = min(R20_low, R20_high);
y_hi = max(R20_low, R20_high);
hR20band = fill(x_band, [y_lo, fliplr(y_hi)], col.r20, ...
                'FaceAlpha',0.20, 'EdgeColor','none');
hR20med  = plot(years_plot, R20_med, 'Color',col.r20, 'LineWidth',2);

% R37 band and medium line
y_lo = min(R37_low, R37_high);
y_hi = max(R37_low, R37_high);
hR37band = fill(x_band, [y_lo, fliplr(y_hi)], col.r37, ...
                'FaceAlpha',0.20, 'EdgeColor','none');
hR37med  = plot(years_plot, R37_med, 'Color',col.r37, 'LineWidth',2);

% R52 band and medium line
y_lo = min(R52_low, R52_high);
y_hi = max(R52_low, R52_high);
hR52band = fill(x_band, [y_lo, fliplr(y_hi)], col.r52, ...
                'FaceAlpha',0.20, 'EdgeColor','none');
hR52med  = plot(years_plot, R52_med, 'Color',col.r52, 'LineWidth',2);

% R100 band and medium line
y_lo = min(R100_low, R100_high);
y_hi = max(R100_low, R100_high);
hR100band = fill(x_band, [y_lo, fliplr(y_hi)], col.r100, ...
                 'FaceAlpha',0.20, 'EdgeColor','none');
hR100med  = plot(years_plot, R100_med, 'Color',col.r100, 'LineWidth',2);

% make sure lines sit above shading
uistack([hR20med hR37med hR52med hR100med],'top');

grid on;
xlabel('Year','FontWeight','bold');
ylabel('Hours of Flooding','FontWeight','bold');
title(sprintf('Annual hours of exceedance %d to %d (%s)', ...
      y_startt, y_endd, loc_title), 'FontWeight','bold');

set(gca,'FontSize',12);

% y axis from zero upward, top auto
ax = gca;
yl = ax.YLim;
ax.YLim = [0 yl(2)];



% x axis limits and ticks tied to user settings
xlim([year_start_plot year_end_plot]);
xticks(year_start_plot:year_tick_step:year_end_plot);
xtickangle(45);
% legend formatted like your example figure
LG = legend([hR20band hR20med ...
        hR37band hR37med ...
        hR52band hR52med ...
        hR100band hR100med], ...
       {'R20 Low-High SLR band','R20 Med SLR', ...
        'R37 Low-High SLR band','R37 Med SLR', ...
        'R52 Low-High SLR band','R52 Med SLR', ...
        'R100 Low-High SLR band','R100 Med SLR'}, ...
       'Location','northwest', 'NumColumns',2);
set(LG,'FontSize',10);   % smaller text
hold off;


% %% ================== Diagnostic M2 and K1 plots for Low Medium High ==================
% % Plot amplitude evolution for M2 and K1 constituents under all SLR scenarios
% %  uses ftidestruc output to display interpolated vs. anchor years
% slr_idx   = [1 2 3];
% slr_tags  = {'l','m','h'};
% slr_names = {'Low SLR','Medium SLR','High SLR'};
% x  = [2020 2040 2060 2080 2100];
% xq = 2020:2100;
% 
% for s = 1:3
%   tag = slr_tags{s};
%   figure('Color','w');
% 
%   %% ---------- M2 ----------
%   subplot(2,1,1)
%   % Plot interpolated amplitude time series for all restoration scenarios
%   plot(xq, squeeze(ftidestruc.r20 .tidecon(11,1,:,slr_idx(s))), 'LineWidth',2, 'Color',col.r20); hold on
%   plot(xq, squeeze(ftidestruc.r37 .tidecon(11,1,:,slr_idx(s))), 'LineWidth',2, 'Color',col.r37);
%   plot(xq, squeeze(ftidestruc.r52 .tidecon(11,1,:,slr_idx(s))), 'LineWidth',2, 'Color',col.r52);
%   plot(xq, squeeze(ftidestruc.r100.tidecon(11,1,:,slr_idx(s))), 'LineWidth',2, 'Color',col.r100);
% 
%   % Anchor points (five-year analyses)
%   plot(x, arrayfun(@(k) T.r20.(tag)(k).tidecon(11,1), 1:5), 'k*','MarkerSize',10)
%   plot(x, arrayfun(@(k) T.r37.(tag)(k).tidecon(11,1), 1:5), 'k*','MarkerSize',10)
%   plot(x, arrayfun(@(k) T.r52.(tag)(k).tidecon(11,1), 1:5), 'k*','MarkerSize',10)
%   plot(x, arrayfun(@(k) T.r100.(tag)(k).tidecon(11,1), 1:5), 'k*','MarkerSize',10)
% 
%   legend('R20','R37','R52','R100','Location','northwest')
%   %title(['M_2  ' slr_names{s}])
%   title (['M_2'])
%   ylim([0.69 0.85])
%   ylabel('Amplitude (m)')
%   set(gca,'FontSize',12)
% 
%   %% ---------- K1 ----------
%   subplot(2,1,2)
%   plot(xq, squeeze(ftidestruc.r20 .tidecon(6,1,:,slr_idx(s))), 'LineWidth',2, 'Color',col.r20); hold on
%   plot(xq, squeeze(ftidestruc.r37 .tidecon(6,1,:,slr_idx(s))), 'LineWidth',2, 'Color',col.r37);
%   plot(xq, squeeze(ftidestruc.r52 .tidecon(6,1,:,slr_idx(s))), 'LineWidth',2, 'Color',col.r52);
%   plot(xq, squeeze(ftidestruc.r100.tidecon(6,1,:,slr_idx(s))), 'LineWidth',2, 'Color',col.r100);
% 
%   plot(x, arrayfun(@(k) T.r20.(tag)(k).tidecon(6,1), 1:5), 'k*','MarkerSize',10)
%   plot(x, arrayfun(@(k) T.r37.(tag)(k).tidecon(6,1), 1:5), 'k*','MarkerSize',10)
%   plot(x, arrayfun(@(k) T.r52.(tag)(k).tidecon(6,1), 1:5), 'k*','MarkerSize',10)
%   plot(x, arrayfun(@(k) T.r100.(tag)(k).tidecon(6,1), 1:5), 'k*','MarkerSize',10)
% 
%   legend('R20','R37','R52','R100','Location','northeast')
%   %title(['K_1  ' slr_names{s}])
%   title (['K_1'])
%   ylim([0.42 0.525])
%   ylabel('Amplitude (m)')
%   set(gca,'FontSize',12)
% end

% ================== Diagnostic M2 K1 and next two tidal constituents for Low Medium High ==================
% M2 and K1 figure kept exactly as before
% Additional figure with two subplots for the next two strongest constituents
% next two are ranked from interpolated ftidestruc amplitudes, excluding M2 and K1

%% ================== Diagnostic M2 K1 and next two tidal constituents for Low Medium High ==================
% idx_M2 = 11;
% idx_K1 = 6;
% 
% amp_all  = squeeze(ftidestruc.r20.tidecon(:,1,:,:));
% amp_mean = mean(amp_all,[2 3]);
% 
% mask = true(size(amp_mean));
% mask([idx_M2 idx_K1]) = false;
% 
% [amp_sorted, order] = sort(amp_mean(mask),'descend');
% all_idx    = find(mask);
% next2_idx  = all_idx(order(1:2));
% name_char  = ftidestruc.r20.name;
% next2_names = cellstr(name_char(next2_idx,:));
% 
% slr_idx   = [1 2 3];
% slr_tags  = {'l','m','h'};
% x  = [2020 2040 2060 2080 2100];
% xq = 2020:2100;
% 
% for s = 1
%   tag = slr_tags{s};
% 
%   % =====================================================
%   % ================== AMPLITUDE FIGURE 1 ================
%   % =====================================================
%   figure('Color','w');
% 
%   %% ---- M2 amplitude ----
%   subplot(2,1,1)
%   plot(xq, squeeze(ftidestruc.r20 .tidecon(11,1,:,slr_idx(s))), 'LineWidth',2, 'Color',col.r20); hold on
%   plot(xq, squeeze(ftidestruc.r37 .tidecon(11,1,:,slr_idx(s))), 'LineWidth',2, 'Color',col.r37);
%   plot(xq, squeeze(ftidestruc.r52 .tidecon(11,1,:,slr_idx(s))), 'LineWidth',2, 'Color',col.r52);
%   plot(xq, squeeze(ftidestruc.r100.tidecon(11,1,:,slr_idx(s))), 'LineWidth',2, 'Color',col.r100);
% 
%   plot(x, arrayfun(@(k) T.r20.(tag)(k).tidecon(11,1), 1:5), 'k*','MarkerSize',10)
%   plot(x, arrayfun(@(k) T.r37.(tag)(k).tidecon(11,1), 1:5), 'k*','MarkerSize',10)
%   plot(x, arrayfun(@(k) T.r52.(tag)(k).tidecon(11,1), 1:5), 'k*','MarkerSize',10)
%   plot(x, arrayfun(@(k) T.r100.(tag)(k).tidecon(11,1), 1:5), 'k*','MarkerSize',10)
% 
%   title('M2')
%   ylim([0.69 0.85])
%   ylabel('Amplitude (m)')
%   set(gca,'FontSize',12)
% 
%   %% ---- K1 amplitude (LEGEND ONLY HERE) ----
%   subplot(2,1,2)
%   plot(xq, squeeze(ftidestruc.r20 .tidecon(6,1,:,slr_idx(s))), 'LineWidth',2, 'Color',col.r20); hold on
%   plot(xq, squeeze(ftidestruc.r37 .tidecon(6,1,:,slr_idx(s))), 'LineWidth',2, 'Color',col.r37);
%   plot(xq, squeeze(ftidestruc.r52 .tidecon(6,1,:,slr_idx(s))), 'LineWidth',2, 'Color',col.r52);
%   plot(xq, squeeze(ftidestruc.r100.tidecon(6,1,:,slr_idx(s))), 'LineWidth',2, 'Color',col.r100);
% 
%   plot(x, arrayfun(@(k) T.r20.(tag)(k).tidecon(6,1), 1:5), 'k*','MarkerSize',10)
%   plot(x, arrayfun(@(k) T.r37.(tag)(k).tidecon(6,1), 1:5), 'k*','MarkerSize',10)
%   plot(x, arrayfun(@(k) T.r52.(tag)(k).tidecon(6,1), 1:5), 'k*','MarkerSize',10)
%   plot(x, arrayfun(@(k) T.r100.(tag)(k).tidecon(6,1), 1:5), 'k*','MarkerSize',10)
% 
%   % ===== dummy star for legend =====
%   h_star = plot(nan, nan, 'k*', 'MarkerSize',10);
% 
%   legend('R20','R37','R52','R100','Hydrodynamic model','Location','northwest')
% 
%   title('K1')
%   ylim([0.42 0.525])
%   ylabel('Amplitude (m)')
%   set(gca,'FontSize',12)
% 
%   % =====================================================
%   % ============ AMPLITUDE FIGURE 2 (next two) ==========
%   % =====================================================
%   figure('Color','w');
% 
%   for p = 1:2
%     ci    = next2_idx(p);
%     cname = strtrim(next2_names{p});
% 
%     subplot(2,1,p); hold on
% 
%     plot(xq, squeeze(ftidestruc.r20 .tidecon(ci,1,:,slr_idx(s))), 'LineWidth',2, 'Color',col.r20);
%     plot(xq, squeeze(ftidestruc.r37 .tidecon(ci,1,:,slr_idx(s))), 'LineWidth',2, 'Color',col.r37);
%     plot(xq, squeeze(ftidestruc.r52 .tidecon(ci,1,:,slr_idx(s))), 'LineWidth',2, 'Color',col.r52);
%     plot(xq, squeeze(ftidestruc.r100.tidecon(ci,1,:,slr_idx(s))), 'LineWidth',2, 'Color',col.r100);
% 
%     plot(x, arrayfun(@(k) T.r20.(tag)(k).tidecon(ci,1), 1:5), 'k*','MarkerSize',10)
%     plot(x, arrayfun(@(k) T.r37.(tag)(k).tidecon(ci,1), 1:5), 'k*','MarkerSize',10)
%     plot(x, arrayfun(@(k) T.r52.(tag)(k).tidecon(ci,1), 1:5), 'k*','MarkerSize',10)
%     plot(x, arrayfun(@(k) T.r100.(tag)(k).tidecon(ci,1), 1:5), 'k*','MarkerSize',10)
%     % 
%     % if p == 2
%     %     h_star = plot(nan, nan, 'k*', 'MarkerSize',10);
%     %     legend('R20','R37','R52','R100','Hydrodynamic model','Location','northwest')
%     % end
% 
%     title(cname)
%     ylabel('Amplitude (m)')
%     set(gca,'FontSize',12)
% 
%     % your double axis overlay (UNCHANGED)
%     ax1 = gca;
%     ax2 = axes('Position',ax1.Position,'Color','none',...
%                'XAxisLocation','top','YAxisLocation','right',...
%                'XLim',ax1.XLim,'YLim',ax1.YLim,...
%                'XTick',ax1.XTick,'YTick',ax1.YTick,...
%                'XTickLabel',[],'YTickLabel',[],...
%                'HitTest','off','HandleVisibility','off');
%     linkaxes([ax1 ax2],'xy');
%   end
% 
%  %=====================================================
%   %%   ===================== PHASE FIGURE 1 ================
%   % =====================================================
%   figure('Color','w');
% 
%   %% ---- M2 phase ----
%   subplot(2,1,1)
%   plot(xq, squeeze(ftidestruc.r20 .tidecon(11,3,:,slr_idx(s))), 'LineWidth',2, 'Color',col.r20); hold on
%   plot(xq, squeeze(ftidestruc.r37 .tidecon(11,3,:,slr_idx(s))), 'LineWidth',2, 'Color',col.r37);
%   plot(xq, squeeze(ftidestruc.r52 .tidecon(11,3,:,slr_idx(s))), 'LineWidth',2, 'Color',col.r52);
%   plot(xq, squeeze(ftidestruc.r100.tidecon(11,3,:,slr_idx(s))), 'LineWidth',2, 'Color',col.r100);
% 
%   plot(x, arrayfun(@(k) T.r20.(tag)(k).tidecon(11,3), 1:5), 'k*','MarkerSize',10)
%   plot(x, arrayfun(@(k) T.r37.(tag)(k).tidecon(11,3), 1:5), 'k*','MarkerSize',10)
%   plot(x, arrayfun(@(k) T.r52.(tag)(k).tidecon(11,3), 1:5), 'k*','MarkerSize',10)
%   plot(x, arrayfun(@(k) T.r100.(tag)(k).tidecon(11,3), 1:5), 'k*','MarkerSize',10)
% 
%   title('M2 phase')
%   ylabel('Phase (deg)')
%   ylim ([200 380])
%   set(gca,'FontSize',12)
% 
%   %% ---- K1 phase (NO LEGEND) ----
%   subplot(2,1,2)
%   plot(xq, squeeze(ftidestruc.r20 .tidecon(6,3,:,slr_idx(s))), 'LineWidth',2, 'Color',col.r20); hold on
%   plot(xq, squeeze(ftidestruc.r37 .tidecon(6,3,:,slr_idx(s))), 'LineWidth',2, 'Color',col.r37);
%   plot(xq, squeeze(ftidestruc.r52 .tidecon(6,3,:,slr_idx(s))), 'LineWidth',2, 'Color',col.r52);
%   plot(xq, squeeze(ftidestruc.r100.tidecon(6,3,:,slr_idx(s))), 'LineWidth',2, 'Color',col.r100);
% 
%   plot(x, arrayfun(@(k) T.r20.(tag)(k).tidecon(6,3), 1:5), 'k*','MarkerSize',10)
%   plot(x, arrayfun(@(k) T.r37.(tag)(k).tidecon(6,3), 1:5), 'k*','MarkerSize',10)
%   plot(x, arrayfun(@(k) T.r52.(tag)(k).tidecon(6,3), 1:5), 'k*','MarkerSize',10)
%   plot(x, arrayfun(@(k) T.r100.(tag)(k).tidecon(6,3), 1:5), 'k*','MarkerSize',10)
% 
%   % *** NO LEGEND ***
%   title('K1 phase')
%   ylabel('Phase (deg)')
%   ylim ([245 280])
%   set(gca,'FontSize',12)
% 
%   % =====================================================
%   % ===================== PHASE FIGURE 2 ================
%   % =====================================================
%   figure('Color','w');
% 
%   for p = 1:2
%     ci    = next2_idx(p);
%     cname = strtrim(next2_names{p});
% 
%     subplot(2,1,p)
%     hold on
% 
%     plot(xq, squeeze(ftidestruc.r20 .tidecon(ci,3,:,slr_idx(s))), 'LineWidth',2, 'Color',col.r20);
%     plot(xq, squeeze(ftidestruc.r37 .tidecon(ci,3,:,slr_idx(s))), 'LineWidth',2, 'Color',col.r37);
%     plot(xq, squeeze(ftidestruc.r52 .tidecon(ci,3,:,slr_idx(s))), 'LineWidth',2, 'Color',col.r52);
%     plot(xq, squeeze(ftidestruc.r100.tidecon(ci,3,:,slr_idx(s))), 'LineWidth',2, 'Color',col.r100);
% 
%     plot(x, arrayfun(@(k) T.r20.(tag)(k).tidecon(ci,3), 1:5), 'k*','MarkerSize',10)
%     plot(x, arrayfun(@(k) T.r37.(tag)(k).tidecon(ci,3), 1:5), 'k*','MarkerSize',10)
%     plot(x, arrayfun(@(k) T.r52.(tag)(k).tidecon(ci,3), 1:5), 'k*','MarkerSize',10)
%     plot(x, arrayfun(@(k) T.r100.(tag)(k).tidecon(ci,3), 1:5), 'k*','MarkerSize',10)
% 
%     % *** NO LEGEND ***
%     title([cname ' phase'])
%     ylabel('Phase (deg)')
%     for p = 1
%         ylim ([200 350]);
%     end
%     set(gca,'FontSize',12)
% 
%     % double axis overlay (UNCHANGED)
%     ax1 = gca;
%     ax2 = axes('Position',ax1.Position,'Color','none',...
%                'XAxisLocation','top','YAxisLocation','right',...
%                'XLim',ax1.XLim,'YLim',ax1.YLim,...
%                'XTick',ax1.XTick,'YTick',ax1.YTick,...
%                'XTickLabel',[],'YTickLabel',[],...
%                'HitTest','off','HandleVisibility','off');
%     linkaxes([ax1 ax2],'xy');
% 
%   end
% 
% end
