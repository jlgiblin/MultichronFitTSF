%% MultichronFitTSF.m
%
% Joint multi-chronometer inversion of detrital age distributions using a
% hypsometry-weighted Gaussian mixture likelihood. All chronometers are
% optimized simultaneously with a closure-temperature-aware ordering penalty
% that enforces the physical age hierarchy (lower-Tc systems yield younger
% ages) without requiring explicit kinetic models.
%
% MODEL DESIGN
%   This script fits all chronometers jointly with fminunc (quasi-Newton),
%   adding a pairwise ordering penalty across all chronometer pairs scaled
%   by their normalized closure temperature separation. This prevents
%   physically implausible age inversions (e.g. ZHe older than APb) while
%   remaining agnostic about the specific thermal history.
%
% USAGE
%   1. Set catchment_name and base_dir below (the only two lines you edit).
%   2. Choose fixed or iterative source weighting in the config file.
%   3. Run. All outputs are written into the catchment subfolder.
%
% FOLDER STRUCTURE
%   <base_dir>/
%   ├── MultichronFitTSF.m
%   ├── SampleA/
%   │   ├── SampleA_config.csv
%   │   ├── SampleA_Hypso.csv
%   │   ├── SampleA_ApHe.csv  (columns: Date_Ma, Error_Ma)
%   │   ├── SampleA_ZHe.csv
%   │   ├── SampleA_ApPb.csv
%   │   ├── SampleA_Hbl.csv
%   │   └── figures_svg/
%   └── ...
%
% CONFIG CSV FORMAT
%   Required columns:
%     Chronometer, File, TauMin, AgeMargin, Lambda, AgeMinFilter,
%     AgeMaxFilter, TSFMode
%   TSFMode must be set to fixed or iterative on the Hypsometry row.
%   Optional global-override columns (add as single-row header + value, or
%   place in the Hypsometry row under these column names):
%     w_order     : ordering penalty weight  (default: see ORDERING PENALTY below)
%     delta_min   : min age separation scale (default: see ORDERING PENALTY below)
%     TSFUpdateFraction : iterative relaxation step [0,1] (default 0.4)
%     TSFSmoothSpan     : iterative smoothing span in bins (default 3)
%   Optional columns on each chronometer row:
%     ClosureTemperature_C : user-defined nominal closure temperature used
%                            only to construct ordering constraints
%     TSFGroup    : matching labels share source weights in iterative mode;
%                   blank means estimate this chronometer separately
%     EstimateTSF : true estimates that group's weights; false keeps the
%                   group fixed to hypsometry (blank defaults to true)
%
%   Example:
%     Chronometer,File,TauMin,AgeMargin,Lambda,AgeMinFilter,AgeMaxFilter,ClosureTemperature_C,TSFMode,TSFGroup,EstimateTSF
%     Hypsometry,SampleA_Hypso.csv,,,,,,,iterative,,
%     ApHe,SampleA_ApHe.csv,0.5,5,0.0,0,Inf,70,,apatite,true
%     ZHe,SampleA_ZHe.csv,0.5,5,0.0,0,Inf,170,,,true
%     ApPb,SampleA_ApPb.csv,0.1,15,0.5,0,Inf,460,,apatite,true
%     Hbl,SampleA_Hbl.csv,0.1,5,1.0,83,Inf,570,,hornblende,false
%
% CLOSURE TEMPERATURES (Tc)
%   Prefer an explicit ClosureTemperature_C value on each chronometer row.
%   Recognized blank entries use built-in fallbacks from Hodges (2014).
%   Unrecognized blank entries are fitted without an ordering constraint.
%
% ORDERING PENALTY
%   For each pair (i,j) with Tc_i < Tc_j, at each elevation bin k:
%     violation = max(0, A_i(k) - A_j(k) + delta_min * gap_ij)^2
%     penalty   = w_order * sum over all pairs and bins of violation
%   where gap_ij = (Tc_j - Tc_i) / (Tc_max - Tc_min).
%   delta_min (Ma) sets the minimum required age separation proportional to
%   the Tc gap -- a proxy for the rapid-cooling end-member.
%   w_order controls penalty stiffness. Both are tunable below.
%
% OUTPUTS  (written to catchment subfolder)
%   predicted_bedrock_transect_<Chron>.csv          (selected solution)
%   predicted_bedrock_transect_<Chron>_CI.csv       (bootstrap CI)
%   source_weighting.csv                            (fixed/iterative weights)
%   source_weighting_convergence.csv                (iterative diagnostics)
%   source_weighting_bootstrap_CI.csv               (iterative uncertainty)
%   source_weighting_bootstrap_summary.csv          (bootstrap diagnostics)
%   grain_posteriors_<Chron>.csv                    (optional)
%   grain_expected_source_<Chron>.csv
%   ordering_violations_preFit.csv                  (pre-fit violations)
%   ordering_penalty_contributions.csv              (per-pair penalty)
%   summary_fit_params.csv
%   figures_svg/
%
% RELEASE
%   v0.1.0: Public fixed/iterative source-weighting modes, explicit optional
%   groups and closure temperatures, reproducible serial/parallel bootstrap
%   uncertainty, source-weight diagnostics, and DEM coordinate assignment.

clear; close all; clc;

%% ============================================================
%% USER SETTINGS  <-- only edit these two lines between runs
% ---- Match these two lines to your MultichronFitTSF.m settings ----
catchment_name = "SampleA";                     % <-- change to your catchment name
base_dir       = "/path/to/MultichronFitTSF";   % <-- change to your project folder
%% ============================================================


%% ---- CLOSURE-TEMPERATURE FALLBACKS ----
%
% Used only when ClosureTemperature_C is blank. Values from Hodges (2014)
% Table 2, bulk closure temperature (Tcb), in degrees Celsius.
%
TC_DEFAULTS = struct( ...
    'ahe',  70, ...   % Apatite (U-Th)/He  -- van Soest et al. 2011; Farley 2000
    'zhe', 170, ...   % Zircon  (U-Th)/He  -- Reiners et al. 2004
    'apb', 460, ...   % Apatite (U-Th)/Pb  -- Cherniak et al. 1991
    'har', 570  ...   % Hornblende 40Ar/39Ar -- Harrison 1981
);
%% -----------------------------------------------------------------


%% ---- ORDERING PENALTY PARAMETERS  <-- tunable ----
%
% w_order    : penalty weight. Higher = stricter age ordering enforcement.
%              Start at 5; increase to 20-50 if violations persist.
%              Can be overridden per-catchment in config CSV (column w_order
%              on the Hypsometry row).
%
% delta_min_Ma : minimum required age separation (Ma) scaled by Tc gap.
%              Represents the rapid-cooling end-member: even in the fastest
%              plausible cooling scenario, ages should differ by at least
%              this amount times the normalised Tc gap.
%              Typical range: 0.5 - 2.0 Ma. Can be overridden in config.
%
w_order      = 5.0;
delta_min_Ma = 1.0;
%% ---------------------------------------------------


%% ---- SOURCE-WEIGHTING MODE ----
%
% TSFMode = "fixed"
%   The measured hypsometry supplies the source weights throughout the fit.
%
% TSFMode = "iterative"
%   Flexible source-weighting analogue inspired by Gallagher & Parra (2020):
%     1. predict a detrital age density for every elevation bin;
%     2. estimate nonnegative, sum-to-one TSF weights by least squares;
%     3. regularize the TSF toward hypsometry and a smooth elevation curve;
%     4. refit A(z)/tau and repeat until the TSF and NLL stabilize.
%   TSFGroup and EstimateTSF in the config file control which chronometers
%   share weights and which groups remain fixed. This estimates effective
%   source contribution by elevation; it does not alter measured hypsometry
%   or reproduce QTQt's thermal-history inversion.
% A mode must be explicitly selected on the Hypsometry row of the config
% file. There is no default source-weighting mode.
%
% tsf_update_fraction controls the departure from hypsometry:
%   0 = fixed hypsometry; 1 = full one-step posterior weight update.
% A value of 0.4 provides a conservative default relaxation step.
% In iterative mode it is the relaxation step from the current weights
% toward the newly estimated regularized nonnegative least-squares target.
%
% tsf_smooth_span is a moving-mean span in equal-area bin index. Use 1 to
% disable smoothing. These values can be overridden on the Hypsometry row
% of the config CSV using TSFMode, TSFUpdateFraction, and TSFSmoothSpan.
tsf_mode            = "";  % required config value; intentionally no default
tsf_update_fraction = 0.4;
tsf_smooth_span     = 3;

% Iterative source-weight settings. These are ignored in fixed mode.
iterative_max_outer       = 40;
iterative_weight_tol      = 1e-4;
iterative_nll_tol         = 1e-3;
iterative_hypsometry_pull = 0.25;
iterative_smoothness      = 1.0;
iterative_no_improve_patience = 3;
iterative_boot_max_outer  = 20;  % per-resample cap (limits bootstrap cost)
%% ----------------------------------------------------------------------


%% ---- Global fitting options (rarely need changing) ----
age_grid_step          = 0.25;   % PDF plot resolution (Ma)
write_grain_posteriors = true;   % set false to skip large grain posterior CSVs
show_optimizer_iters   = false;  % show fminunc iterations in console
do_bootstrap           = false;  % enable after checking the primary fit
n_boot                 = 20;     % pilot resamples (about 200 for a final run)
bootstrap_random_seed  = 1;      % nonnegative integer for reproducible resampling
use_parallel_bootstrap = false;  % true uses Parallel Computing Toolbox when available
parallel_worker_count  = 4;      % [] uses the local profile default
ci_lo                  = 0.16;   % 68% CI lower bound
ci_hi                  = 0.84;   % 68% CI upper bound
target_bins            = 20;     % equal-area hypsometry bins ([] = use raw)
min_bins               = 8;      % safety floor on bin count


%% ---------------- LOCATE AND READ CONFIG FILE ----------------
catchment_dir = fullfile(base_dir, catchment_name);
config_file   = fullfile(catchment_dir, catchment_name + "_config.csv");

if ~isfile(config_file)
    error("Config file not found: %s\n" + ...
          "Expected: <base_dir>/<catchment_name>/<catchment_name>_config.csv", ...
          config_file);
end

cfg = readtable(config_file, 'TextType', 'string');

% Validate required columns
required_cols = {'Chronometer','File','TauMin','AgeMargin', ...
                 'Lambda','AgeMinFilter','AgeMaxFilter','TSFMode'};
for rc = required_cols
    if ~any(strcmpi(cfg.Properties.VariableNames, rc{1}))
        error("Config file is missing required column: %s", rc{1});
    end
end

fprintf("Loaded config : %s\n",   config_file);
fprintf("Catchment     : %s\n",   catchment_name);
fprintf("Catchment dir : %s\n\n", catchment_dir);

% ---- Read global settings from the Hypsometry row ----
hyps_mask  = strcmpi(cfg.Chronometer, 'Hypsometry');
chron_mask = ~hyps_mask;

if sum(hyps_mask) ~= 1
    error("Config must contain exactly one row with Chronometer = 'Hypsometry'.");
end

hyps_row = cfg(hyps_mask, :);
if any(strcmpi(cfg.Properties.VariableNames, 'w_order')) && ...
        ~isnan(double(hyps_row.w_order))
    w_order = double(hyps_row.w_order);
    fprintf("Config override: w_order = %.2f\n", w_order);
end
if any(strcmpi(cfg.Properties.VariableNames, 'delta_min')) && ...
        ~isnan(double(hyps_row.delta_min))
    delta_min_Ma = double(hyps_row.delta_min);
    fprintf("Config override: delta_min_Ma = %.2f\n", delta_min_Ma);
end

tsf_mode_cfg = string(hyps_row.TSFMode);
if ismissing(tsf_mode_cfg) || strlength(strtrim(tsf_mode_cfg)) == 0
    error("TSFMode must be explicitly set to 'fixed' or 'iterative' " + ...
          "on the Hypsometry row. There is no default mode.");
end
tsf_mode = lower(strtrim(tsf_mode_cfg));
fprintf("Config setting: TSFMode = %s\n", tsf_mode);
if any(strcmpi(cfg.Properties.VariableNames, 'TSFUpdateFraction'))
    tsf_alpha_cfg = double(hyps_row.TSFUpdateFraction);
    if ~isnan(tsf_alpha_cfg)
        tsf_update_fraction = tsf_alpha_cfg;
        fprintf("Config override: TSFUpdateFraction = %.3f\n", ...
            tsf_update_fraction);
    end
end
if any(strcmpi(cfg.Properties.VariableNames, 'TSFSmoothSpan'))
    tsf_span_cfg = double(hyps_row.TSFSmoothSpan);
    if ~isnan(tsf_span_cfg)
        tsf_smooth_span = tsf_span_cfg;
        fprintf("Config override: TSFSmoothSpan = %d\n", tsf_smooth_span);
    end
end

if ~any(tsf_mode == ["fixed", "iterative"])
    error("TSFMode must be 'fixed' or 'iterative'; received '%s'.", tsf_mode);
end
if ~isscalar(tsf_update_fraction) || ~isfinite(tsf_update_fraction) || ...
        tsf_update_fraction < 0 || tsf_update_fraction > 1
    error("TSFUpdateFraction must be a finite scalar between 0 and 1.");
end
if ~isscalar(tsf_smooth_span) || ~isfinite(tsf_smooth_span) || ...
        tsf_smooth_span < 1 || tsf_smooth_span ~= round(tsf_smooth_span)
    error("TSFSmoothSpan must be a positive integer.");
end
tsf_smooth_span = round(tsf_smooth_span);

if ~isscalar(iterative_max_outer) || iterative_max_outer < 1 || ...
        iterative_max_outer ~= round(iterative_max_outer)
    error("iterative_max_outer must be a positive integer.");
end
if ~isscalar(iterative_no_improve_patience) || ...
        iterative_no_improve_patience < 1 || ...
        iterative_no_improve_patience ~= round(iterative_no_improve_patience)
    error("iterative_no_improve_patience must be a positive integer.");
end
if ~isscalar(iterative_boot_max_outer) || iterative_boot_max_outer < 1 || ...
        iterative_boot_max_outer ~= round(iterative_boot_max_outer)
    error("iterative_boot_max_outer must be a positive integer.");
end
if any(~isfinite([iterative_weight_tol, iterative_nll_tol, ...
                  iterative_hypsometry_pull, iterative_smoothness])) || ...
        any([iterative_weight_tol, iterative_nll_tol, ...
             iterative_hypsometry_pull, iterative_smoothness] < 0)
    error("Iterative source-weight tolerances and regularization settings must be finite and nonnegative.");
end
if ~isscalar(n_boot) || ~isfinite(n_boot) || n_boot < 1 || ...
        n_boot ~= round(n_boot)
    error("n_boot must be a positive integer.");
end
if ~isscalar(bootstrap_random_seed) || ~isfinite(bootstrap_random_seed) || ...
        bootstrap_random_seed < 0 || ...
        bootstrap_random_seed > 2^32 - 1 || ...
        bootstrap_random_seed ~= round(bootstrap_random_seed)
    error("bootstrap_random_seed must be an integer from 0 through 2^32-1.");
end
if ~(islogical(use_parallel_bootstrap) || isnumeric(use_parallel_bootstrap)) || ...
        ~isscalar(use_parallel_bootstrap)
    error("use_parallel_bootstrap must be true or false.");
end
use_parallel_bootstrap = logical(use_parallel_bootstrap);
if ~isempty(parallel_worker_count) && ...
        (~isscalar(parallel_worker_count) || ~isfinite(parallel_worker_count) || ...
         parallel_worker_count < 1 || parallel_worker_count ~= round(parallel_worker_count))
    error("parallel_worker_count must be empty or a positive integer.");
end
hypsometry_file = fullfile(catchment_dir, cfg.File(hyps_mask));
chron_cfg       = cfg(chron_mask, :);
n_chron         = height(chron_cfg);

fprintf("Hypsometry file : %s\n", hypsometry_file);
fprintf("Chronometers    : %s\n", strjoin(chron_cfg.Chronometer, ', '));
fprintf("w_order         : %.2f\n", w_order);
fprintf("delta_min_Ma    : %.2f\n", delta_min_Ma);
fprintf("TSF mode        : %s\n", tsf_mode);
if tsf_mode == "iterative"
    fprintf("Iterative update : relaxation=%.2f  smooth span=%d bins\n", ...
        tsf_update_fraction, tsf_smooth_span);
    fprintf("Regularization   : hypsometry=%.3g  smoothness=%.3g  max outer=%d\n", ...
        iterative_hypsometry_pull, iterative_smoothness, iterative_max_outer);
    fprintf("Stopping rule    : stop after %d consecutive non-improving iterations\n", ...
        iterative_no_improve_patience);
    if do_bootstrap
        fprintf("Iterative bootstrap: re-estimate source weights, max outer=%d per resample\n", ...
            iterative_boot_max_outer);
    end
end
fprintf("\n");


%% ---------------- READ AND PROCESS HYPSOMETRY ----------------
if ~isfile(hypsometry_file)
    error("Hypsometry file not found: %s", hypsometry_file);
end

H = readtable(hypsometry_file);

% Detect elevation column
if any(strcmpi(H.Properties.VariableNames, 'Elevation_m'))
    z_raw = H.Elevation_m(:);
elseif any(strcmpi(H.Properties.VariableNames, 'Elevation'))
    z_raw = H.Elevation(:);
else
    error("Hypsometry file must contain 'Elevation_m' or 'Elevation'.");
end

% Detect CDF column
hasRel  = any(strcmpi(H.Properties.VariableNames, 'RelArea'));
hasArea = any(strcmpi(H.Properties.VariableNames, 'Area'));
if hasRel
    F_raw = H.RelArea(:);
elseif hasArea
    A_raw = H.Area(:);
    if any(diff(A_raw) < -1e-10)
        error("Hypsometry 'Area' must be nondecreasing with elevation.");
    end
    if A_raw(end) <= 0; error("Hypsometry 'Area' final value must be > 0."); end
    F_raw = A_raw / A_raw(end);
else
    error("Hypsometry file must contain 'RelArea' or 'Area'.");
end

% Sort, deduplicate, validate
[z_raw, sort_idx] = sort(z_raw);
F_raw = F_raw(sort_idx);
[zu, ia] = unique(z_raw, 'stable');
Fu = F_raw(ia);
if any(diff(Fu) < -1e-10); error('Hypsometry CDF must be nondecreasing.'); end
if Fu(end) == Fu(1);       error('Hypsometry CDF is constant.'); end

% Normalise to [0,1]
Fu = (Fu - Fu(1)) / (Fu(end) - Fu(1));
Fu(1) = 0; Fu(end) = 1;

% Resample to equal-area bins
if ~isempty(target_bins)
    target_bins = max(target_bins, min_bins);
    Fq = linspace(0, 1, target_bins + 1)';
    [Fu_u, iu] = unique(Fu, 'stable');
    zu_u = zu(iu);
    z_edges = interp1(Fu_u, zu_u, Fq, 'linear', 'extrap');
    for k = 2:numel(z_edges)
        z_edges(k) = max(z_edges(k), z_edges(k-1));
    end
    z_edges(1) = zu_u(1); z_edges(end) = zu_u(end);
    relArea = Fq;
else
    z_edges = zu;
    relArea = Fu;
end

pz = diff(relArea);
pz(pz < 0) = 0;
pz = pz / sum(pz);

z_centers = 0.5 * (z_edges(1:end-1) + z_edges(2:end));
Nz = numel(z_centers);

fprintf("Hypsometry: %d equal-area bins  |  %.0f - %.0f m\n\n", ...
    Nz, z_edges(1), z_edges(end));


%% ---------------- BUILD Tc LOOKUP AND ASSIGN Tc PER CHRONOMETER ----
tc_lookup = build_tc_lookup(TC_DEFAULTS);

Tc_vec = nan(n_chron, 1);   % closure temperature for each chronometer
Tc_source = repmat("not provided", n_chron, 1);
tc_col = find(strcmpi(chron_cfg.Properties.VariableNames, ...
    'ClosureTemperature_C'), 1);
for c = 1:n_chron
    if ~isempty(tc_col)
        tc_config_text = strtrim(string(chron_cfg{c, tc_col}));
        if ~ismissing(tc_config_text) && strlength(tc_config_text) > 0
            tc_config_value = str2double(tc_config_text);
            if ~isfinite(tc_config_value)
                error("ClosureTemperature_C for chronometer '%s' must be a finite number or blank.", ...
                    chron_cfg.Chronometer(c));
            end
            Tc_vec(c) = tc_config_value;
            Tc_source(c) = "config";
        end
    end

    if isfinite(Tc_vec(c))
        continue;
    end

    cname = lower(strtrim(char(chron_cfg.Chronometer(c))));
    for row = 1:size(tc_lookup, 1)
        if any(strcmpi(tc_lookup{row,1}, cname))
            Tc_vec(c) = tc_lookup{row,2};
            Tc_source(c) = "built-in fallback";
            break;
        end
    end
    if isnan(Tc_vec(c))
        warning("No Tc match for chronometer '%s'. " + ...
            "It will be excluded from ordering penalty.", ...
            chron_cfg.Chronometer(c));
    end
end

% Chronometers with a recognised Tc participate in ordering penalty
has_tc = ~isnan(Tc_vec);
fprintf("Tc assignments:\n");
for c = 1:n_chron
    if has_tc(c)
        fprintf("  %-12s  Tc = %g degC (%s)\n", ...
            chron_cfg.Chronometer(c), Tc_vec(c), Tc_source(c));
    else
        fprintf("  %-12s  Tc = unrecognised (no ordering constraint)\n", ...
            chron_cfg.Chronometer(c));
    end
end

% Build pairwise gap matrix: gap(i,j) = (Tc_j - Tc_i) / (Tc_max - Tc_min)
% Only populated where i < j and Tc_i < Tc_j (lower-Tc must be younger).
Tc_known  = Tc_vec(has_tc);
if numel(Tc_known) >= 2
    Tc_range = max(Tc_known) - min(Tc_known);
else
    Tc_range = 1;
end
if Tc_range == 0; Tc_range = 1; end

% Full index pairs (using original chronometer indices)
chron_idx_with_tc = find(has_tc);
n_pairs = 0;
pair_lo = [];  % index of lower-Tc chronometer in each pair
pair_hi = [];  % index of higher-Tc chronometer
pair_gap = []; % normalised Tc gap [0,1]

for ii = 1:numel(chron_idx_with_tc)
    for jj = ii+1:numel(chron_idx_with_tc)
        ci = chron_idx_with_tc(ii);
        cj = chron_idx_with_tc(jj);
        if Tc_vec(ci) == Tc_vec(cj)
            continue;
        elseif Tc_vec(ci) < Tc_vec(cj)
            lo = ci; hi = cj;
        else
            lo = cj; hi = ci;
        end
        n_pairs  = n_pairs + 1;
        pair_lo(end+1) = lo;  %#ok<AGROW>
        pair_hi(end+1) = hi;  %#ok<AGROW>
        pair_gap(end+1) = abs(Tc_vec(ci) - Tc_vec(cj)) / Tc_range; %#ok<AGROW>
    end
end

fprintf("\nOrdering pairs (%d total):\n", n_pairs);
for p = 1:n_pairs
    fprintf("  %s (Tc=%g) < %s (Tc=%g)  |  gap=%.3f  |  min_sep=%.3f Ma\n", ...
        chron_cfg.Chronometer(pair_lo(p)), Tc_vec(pair_lo(p)), ...
        chron_cfg.Chronometer(pair_hi(p)), Tc_vec(pair_hi(p)), ...
        pair_gap(p), delta_min_Ma * pair_gap(p));
end
fprintf("\n");


%% ---------------- READ ALL DETRITAL DATA ----------------
% Load all chronometer data first so the joint objective can access them.
fprintf("Loading detrital data...\n");
chron_data = struct();   % chron_data(c).age, .sig, .N, .AminB, .AmaxB, etc.
tsf_group_col = find(strcmpi(chron_cfg.Properties.VariableNames, 'TSFGroup'), 1);
estimate_tsf_col = find(strcmpi(chron_cfg.Properties.VariableNames, 'EstimateTSF'), 1);

for c = 1:n_chron
    chron     = char(chron_cfg.Chronometer(c));
    data_file = fullfile(catchment_dir, chron_cfg.File(c));

    chron_data(c).chron        = chron;
    chron_data(c).tau_min      = chron_cfg.TauMin(c);
    chron_data(c).age_margin   = chron_cfg.AgeMargin(c);
    chron_data(c).lambda       = chron_cfg.Lambda(c);
    chron_data(c).age_min_filt = chron_cfg.AgeMinFilter(c);
    chron_data(c).age_max_filt = chron_cfg.AgeMaxFilter(c);
    if isnan(chron_data(c).age_max_filt) || isinf(chron_data(c).age_max_filt)
        chron_data(c).age_max_filt = Inf;
    end
    if tsf_mode == "iterative"
        % Blank or absent TSFGroup means this chronometer gets its own
        % independent source-weight curve. Matching labels opt into sharing.
        chron_data(c).tsf_group = string(chron);
        if ~isempty(tsf_group_col)
            configured_group = strtrim(string(chron_cfg{c, tsf_group_col}));
            if ~ismissing(configured_group) && strlength(configured_group) > 0
                chron_data(c).tsf_group = configured_group;
            end
        end

        % Iterative estimation is the default. Users can explicitly hold a
        % chronometer or group fixed to hypsometry with EstimateTSF=false.
        chron_data(c).estimate_tsf = true;
        if ~isempty(estimate_tsf_col)
            configured_estimate = string(chron_cfg{c, estimate_tsf_col});
            if ~ismissing(configured_estimate) && ...
                    strlength(strtrim(configured_estimate)) > 0
                chron_data(c).estimate_tsf = parse_config_logical( ...
                    chron_cfg{c, estimate_tsf_col}, "EstimateTSF", chron);
            end
        end
    else
        chron_data(c).tsf_group = "fixed_hypsometry";
        chron_data(c).estimate_tsf = false;
    end
    chron_data(c).ok = false;

    if ~isfile(data_file)
        warning("Data file not found: %s. Skipping %s.", data_file, chron);
        continue;
    end

    D   = readtable(data_file);
    age = D.Date_Ma(:);
    sig = D.Error_Ma(:);

    mask = ~isnan(age) & ~isnan(sig) & sig > 0;
    age  = age(mask);
    sig  = sig(mask);

    filt     = age >= chron_data(c).age_min_filt & ...
               age <= chron_data(c).age_max_filt;
    n_excl   = sum(~filt);
    age      = age(filt);
    sig      = sig(filt);
    N        = numel(age);

    if N < 10
        warning("Too few grains (%d) for %s after filtering. Skipping.", N, chron);
        continue;
    end

    chron_data(c).age        = age;
    chron_data(c).sig        = sig;
    chron_data(c).N          = N;
    chron_data(c).n_excluded = n_excl;
    chron_data(c).AminB      = max(0, min(age) - chron_data(c).age_margin);
    chron_data(c).AmaxB      = max(age) + chron_data(c).age_margin;
    chron_data(c).ok         = true;

    fprintf("  %-12s  N=%d  excluded=%d  age=[%.1f, %.1f] Ma\n", ...
        chron, N, n_excl, min(age), max(age));
end
fprintf("\n");

% Identify which chronometers loaded successfully
ok_mask = arrayfun(@(s) s.ok, chron_data);
ok_idx  = find(ok_mask);
n_ok    = numel(ok_idx);

if n_ok == 0
    error("No chronometers loaded successfully. Check data files.");
end


%% ---------------- INDEPENDENT INITIALISATION ----------------
% Fit each chronometer independently (fminsearch, fast) to get good
% starting parameters for the joint fminunc solve.

fprintf("Phase 1: independent initialisation (fminsearch)...\n");

opts_init          = optimset('MaxIter', 1e4, 'MaxFunEvals', 1e5, 'Display', 'off');
theta_init_all     = cell(n_chron, 1);  % starting theta per chronometer
Ahat_init          = nan(Nz, n_chron);  % A(z) from independent fits
tauhat_init        = nan(1,  n_chron);

for c = ok_idx
    cd   = chron_data(c);
    age  = cd.age;  sig = cd.sig;
    AminB = cd.AminB;  AmaxB = cd.AmaxB;
    tau_min = cd.tau_min;

    % Monotonic initial guess via hypsometry-CDF mapping
    Fz = cumsum(pz);
    as = sort(age);
    Fa = ((1:numel(as))' - 0.5) / numel(as);
    A0 = interp1(Fa, as, Fz, 'linear', 'extrap');
    A0 = max(AminB, min(AmaxB, A0));
    for k = 2:Nz; A0(k) = max(A0(k), A0(k-1)); end

    epsv = 1e-6;
    q    = (A0 - AminB) / (AmaxB - AminB);
    q    = max(epsv, min(1-epsv, q));
    s0   = log(q ./ (1-q));
    tau0 = max(tau_min + 0.1, median(sig));
    theta0 = [s0; log(tau0 - tau_min)];

    nll_fn  = @(th) nll_single(th, age, sig, pz, AminB, AmaxB, tau_min, cd.lambda);
    th_hat  = fminsearch(nll_fn, theta0, opts_init);

    [Ahat_init(:,c), tauhat_init(c)] = unpack_single(th_hat, Nz, AminB, AmaxB, tau_min);
    theta_init_all{c} = th_hat;

    fprintf("  %-12s  init tau=%.3f  A=[%.1f, %.1f] Ma\n", ...
        cd.chron, tauhat_init(c), min(Ahat_init(:,c)), max(Ahat_init(:,c)));
end
fprintf("\n");


%% ---------------- REPORT PRE-FIT ORDERING VIOLATIONS ----------------
fprintf("Pre-fit ordering violations:\n");
ViolTable = table();
any_viol  = false;
for p = 1:n_pairs
    lo = pair_lo(p);  hi = pair_hi(p);
    if ~(ok_mask(lo) && ok_mask(hi)); continue; end
    diff_vec   = Ahat_init(:,lo) - Ahat_init(:,hi);  % should be negative
    min_sep    = delta_min_Ma * pair_gap(p);
    viol_mask  = diff_vec > -min_sep;
    if any(viol_mask)
        any_viol = true;
        for k = find(viol_mask)'
            fprintf("  elev=%.0fm  %s(%.2f) > %s(%.2f)  excess=%.3f Ma\n", ...
                z_centers(k), ...
                chron_data(lo).chron, Ahat_init(k,lo), ...
                chron_data(hi).chron, Ahat_init(k,hi), ...
                diff_vec(k) + min_sep);
            row = table(z_centers(k), ...
                string(chron_data(lo).chron), Ahat_init(k,lo), ...
                string(chron_data(hi).chron), Ahat_init(k,hi), ...
                diff_vec(k) + min_sep, ...
                'VariableNames', {'Elevation_m','Chron_Lo','Age_Lo_Ma', ...
                                  'Chron_Hi','Age_Hi_Ma','Violation_Ma'});
            ViolTable = [ViolTable; row]; %#ok<AGROW>
        end
    end
end
if ~any_viol
    fprintf("  None detected (independent fits already satisfy ordering).\n");
end
if ~isempty(ViolTable)
    viol_file = fullfile(catchment_dir, "ordering_violations_preFit.csv");
    writetable(ViolTable, viol_file);
    fprintf("  Wrote %s\n", viol_file);
end
fprintf("\n");


%% ---------------- BUILD JOINT THETA0 ----------------
% Concatenate [s_c(1..Nz); logtau_c] for each ok chronometer, in order.
% Chronometers that failed to load are excluded from the joint vector.

theta0_joint = [];
for c = ok_idx
    theta0_joint = [theta0_joint; theta_init_all{c}]; %#ok<AGROW>
end
% theta0_joint has length n_ok * (Nz+1)


%% ---------------- JOINT fminunc OPTIMISATION ----------------
fprintf("Phase 2: joint fminunc optimisation (%d chronometers, %d params)...\n", ...
    n_ok, numel(theta0_joint));

if show_optimizer_iters
    disp_opt = 'iter';
else
    disp_opt = 'off';
end

opts_joint = optimoptions('fminunc', ...
    'Algorithm',           'quasi-newton', ...
    'MaxIterations',       5000, ...
    'MaxFunctionEvaluations', 1e6, ...
    'FiniteDifferenceType','central', ...
    'Display',             disp_opt, ...
    'OptimalityTolerance', 1e-8, ...
    'StepTolerance',       1e-10);

joint_obj = @(th) nll_joint(th, chron_data, ok_idx, pz, Nz, ...
                             pair_lo, pair_hi, pair_gap, ...
                             w_order, delta_min_Ma);

[theta_joint_hat, nll_total] = fminunc(joint_obj, theta0_joint, opts_joint);

fprintf("Fixed-hypsometry joint fit complete.  Total NLL = %.4f\n\n", nll_total);

% Preserve the fixed-hypsometry reference solution. Iterative-mode outputs
% retain it for direct comparison with the selected source-weight solution.
theta_fixed_hat = theta_joint_hat;
nll_fixed       = nll_total;
[Ahat_all, tauhat_all] = unpack_joint(theta_joint_hat, chron_data, ok_idx, Nz);
Ahat_fixed   = Ahat_all;
tauhat_fixed = tauhat_all;


%% ---------------- OPTIONAL TSF REWEIGHTING ----------------
% The raw update is retained as a diagnostic. Iterative mode alternates a
% regularized nonnegative least-squares source-weight update with A(z)/tau
% refits. Fixed mode retains measured hypsometric weights throughout.

[tsf_raw_update, n_tsf_grains] = posterior_tsf_update( ...
    Ahat_fixed, tauhat_fixed, chron_data, ok_idx, pz);

tsf_weights = pz;
source_weight_raw = nan(Nz, 1);
chron_tsf_group = repmat("fixed_hypsometry", n_chron, 1);
TSFConvergence = table();
% Bootstrap helper arguments are unused in fixed mode, but initialize them
% so fixed-source-weight runs can enter the shared bootstrap code path.
tsf_group_indices = {};
tsf_group_fixed = false(0, 1);

if tsf_mode == "iterative"
    fprintf("Phase 2b: iterative source-weight inversion...\n");

    [tsf_group_names, tsf_group_indices, tsf_group_fixed, ...
        chron_tsf_group] = build_configured_tsf_groups(chron_data, ok_idx);
    weights_current = repmat(pz, 1, n_chron);
    raw_current = repmat(pz, 1, n_chron);
    fprintf("  Source-weight groups:\n");
    for g = 1:numel(tsf_group_names)
        members = tsf_group_indices{g};
        member_names = string(arrayfun(@(c) chron_data(c).chron, ...
            members, 'UniformOutput', false));
        fprintf("    %-18s : %s", tsf_group_names(g), ...
            strjoin(member_names, ', '));
        if tsf_group_fixed(g); fprintf(" (fixed to hypsometry)"); end
        fprintf("\n");
    end

    theta_current = theta_fixed_hat;
    nll_previous = nll_fixed;

    % Retain the best-likelihood state. The source-weight density objective and the
    % grain likelihood are related but not identical, so alternating updates
    % are not guaranteed to improve both.
    best_theta       = theta_fixed_hat;
    best_weights     = weights_current;
    best_source_weight_raw = raw_current;
    best_nll         = nll_fixed;
    best_iteration   = 0;
    no_improve_count = 0;
    stopped_no_improve = false;

    hist_iter       = nan(iterative_max_outer, 1);
    hist_nll        = nan(iterative_max_outer, 1);
    hist_nll_change = nan(iterative_max_outer, 1);
    hist_max_dw     = nan(iterative_max_outer, 1);
    hist_tv         = nan(iterative_max_outer, 1);
    hist_weight_sse  = nan(iterative_max_outer, 1);
    hist_weight_rows = nan(iterative_max_outer, 1);
    converged = false;

    for outer = 1:iterative_max_outer
        [A_current, tau_current] = unpack_joint( ...
            theta_current, chron_data, ok_idx, Nz);

        weights_new = weights_current;
        raw_new     = raw_current;
        weight_fit_sse  = 0;
        weight_fit_rows = 0;

        for g = 1:numel(tsf_group_names)
            members = tsf_group_indices{g};
            if tsf_group_fixed(g)
                raw_g = pz;
                target_g = pz;
                sse_g = 0;
                rows_g = 0;
            else
                [raw_g, sse_g, rows_g] = estimate_source_weights_nnls( ...
                    A_current, tau_current, chron_data, members, pz, ...
                    age_grid_step, iterative_hypsometry_pull, ...
                    iterative_smoothness);
                target_g = raw_g;
                if tsf_smooth_span > 1
                    target_g = smoothdata(target_g, 'movmean', tsf_smooth_span);
                end
                target_g = max(target_g, 1e-8);
                target_g = target_g / sum(target_g);
            end

            for c_group = members
                raw_new(:,c_group) = raw_g;
                if tsf_group_fixed(g)
                    weights_new(:,c_group) = pz;
                else
                    weights_new(:,c_group) = ...
                        (1 - tsf_update_fraction) * weights_current(:,c_group) + ...
                         tsf_update_fraction * target_g;
                    weights_new(:,c_group) = max(weights_new(:,c_group), 1e-8);
                    weights_new(:,c_group) = weights_new(:,c_group) / ...
                        sum(weights_new(:,c_group));
                end
            end
            weight_fit_sse  = weight_fit_sse + sse_g;
            weight_fit_rows = weight_fit_rows + rows_g;
        end

        iterative_obj = @(th) nll_joint(th, chron_data, ok_idx, weights_new, Nz, ...
                                        pair_lo, pair_hi, pair_gap, ...
                                        w_order, delta_min_Ma);
        [theta_new, nll_new] = fminunc(iterative_obj, theta_current, opts_joint);

        max_dw = max(abs(weights_new(:) - weights_current(:)));
        nll_change = nll_new - nll_previous;
        tv_by_group = zeros(numel(tsf_group_names), 1);
        for g = 1:numel(tsf_group_names)
            member0 = tsf_group_indices{g}(1);
            tv_by_group(g) = 0.5 * sum(abs(weights_new(:,member0) - pz));
        end
        tv_hypso = max(tv_by_group);

        hist_iter(outer)       = outer;
        hist_nll(outer)        = nll_new;
        hist_nll_change(outer) = nll_change;
        hist_max_dw(outer)     = max_dw;
        hist_tv(outer)         = tv_hypso;
        hist_weight_sse(outer)  = weight_fit_sse;
        hist_weight_rows(outer) = weight_fit_rows;

        fprintf("  outer %2d: NLL=%.4f  dNLL=%+.4g  max|dw|=%.3g  TV=%.4f\n", ...
            outer, nll_new, nll_change, max_dw, tv_hypso);

        if nll_new < best_nll
            best_theta       = theta_new;
            best_weights     = weights_new;
            best_source_weight_raw = raw_new;
            best_nll         = nll_new;
            best_iteration   = outer;
            no_improve_count = 0;
        else
            no_improve_count = no_improve_count + 1;
        end

        theta_current   = theta_new;
        weights_current = weights_new;
        raw_current     = raw_new;
        nll_previous    = nll_new;

        if outer >= 2 && max_dw <= iterative_weight_tol && ...
                abs(nll_change) <= iterative_nll_tol
            converged = true;
            break;
        end
        if no_improve_count >= iterative_no_improve_patience
            stopped_no_improve = true;
            fprintf("  stopping: %d consecutive iterations failed to improve NLL.\n", ...
                iterative_no_improve_patience);
            break;
        end
    end

    n_outer = find(~isnan(hist_iter), 1, 'last');
    theta_joint_hat = best_theta;
    tsf_weights  = best_weights;
    source_weight_raw = best_source_weight_raw;
    nll_total       = best_nll;
    [Ahat_all, tauhat_all] = unpack_joint( ...
        theta_joint_hat, chron_data, ok_idx, Nz);

    TSFConvergence = table( ...
        hist_iter(1:n_outer), hist_nll(1:n_outer), ...
        hist_nll_change(1:n_outer), hist_max_dw(1:n_outer), ...
        hist_tv(1:n_outer), hist_weight_sse(1:n_outer), ...
        hist_weight_rows(1:n_outer), ...
        'VariableNames', {'Iteration','TotalNLL','NLLChange', ...
                          'MaxAbsoluteWeightChange', ...
                          'TotalVariationFromHypsometry', ...
                          'SourceWeightDataSSE','SourceWeightDataRows'});
    TSFConvergence.Converged = repmat(converged, n_outer, 1);
    TSFConvergence.IsSelectedBest = hist_iter(1:n_outer) == best_iteration;
    if converged
        termination_reason = "converged_tolerances";
    elseif stopped_no_improve
        termination_reason = "stopped_no_likelihood_improvement";
    else
        termination_reason = "maximum_outer_iterations";
    end
    TSFConvergence.TerminationReason = repmat(termination_reason, n_outer, 1);

    conv_file = fullfile(catchment_dir, "source_weighting_convergence.csv");
    writetable(TSFConvergence, conv_file);
    fprintf("Iterative source-weight fit complete after %d evaluated outer iterations.\n", n_outer);
    fprintf("  Selected iteration=%d  termination=%s\n", ...
        best_iteration, termination_reason);
    fprintf("  Total NLL=%.4f (change from fixed=%+.4f)\n", ...
        nll_total, nll_total - nll_fixed);
    for g = 1:numel(tsf_group_names)
        member0 = tsf_group_indices{g}(1);
        group_tv = 0.5 * sum(abs(tsf_weights(:,member0) - pz));
        fprintf("  Source-weight group %-18s TV from hypsometry=%.4f\n", ...
            tsf_group_names(g), group_tv);
    end
    fprintf("  Wrote %s\n\n", conv_file);
else
    fprintf("TSF remains fixed to hypsometry. Raw posterior update retained " + ...
            "for diagnostics only.\n\n");
end

if tsf_mode == "iterative"
    TSFTable = table();
    for g = 1:numel(tsf_group_names)
        members = tsf_group_indices{g};
        member0 = members(1);
        w_g = tsf_weights(:,member0);
        raw_g = source_weight_raw(:,member0);
        [post_g, n_group_grains] = posterior_tsf_update( ...
            Ahat_fixed, tauhat_fixed, chron_data, members, pz);
        group_col = repmat(tsf_group_names(g), Nz, 1);
        fixed_col = repmat(tsf_group_fixed(g), Nz, 1);
        member_col = repmat(strjoin(string(arrayfun( ...
            @(c) chron_data(c).chron, members, 'UniformOutput', false)), '+'), Nz, 1);
        row_group = table(group_col, member_col, fixed_col, z_centers, pz, ...
            post_g, raw_g, w_g, w_g ./ max(pz, realmin), ...
            cumsum(pz), cumsum(w_g), repmat(n_group_grains, Nz, 1), ...
            'VariableNames', {'TSFGroup','Chronometers','FixedToHypsometry', ...
                              'Elevation_m','HypsometryWeight', ...
                              'RawPosteriorWeightUpdate','RawEstimatedWeightUpdate', ...
                              'ModelTSFWeight','RelativeYieldVsHypsometry', ...
                              'HypsometryCDF','ModelTSFCDF','GroupGrainCount'});
        TSFTable = [TSFTable; row_group]; %#ok<AGROW>
    end
    tsf_file = fullfile(catchment_dir, "source_weighting.csv");
else
    tsf_relative_yield = tsf_weights ./ max(pz, realmin);
    TSFTable = table(z_centers, pz, tsf_raw_update, source_weight_raw, tsf_weights, ...
        tsf_relative_yield, cumsum(pz), cumsum(tsf_weights), ...
        'VariableNames', {'Elevation_m','HypsometryWeight', ...
                          'RawPosteriorWeightUpdate','RawEstimatedWeightUpdate', ...
                          'ModelTSFWeight', ...
                          'RelativeYieldVsHypsometry','HypsometryCDF','ModelTSFCDF'});
    tsf_file = fullfile(catchment_dir, "source_weighting.csv");
end

TSFTable.TSFMode = repmat(tsf_mode, height(TSFTable), 1);
TSFTable.UpdateFraction = repmat(tsf_update_fraction, height(TSFTable), 1);
TSFTable.SmoothSpan = repmat(tsf_smooth_span, height(TSFTable), 1);
TSFTable.HypsometryPull = repmat(iterative_hypsometry_pull, height(TSFTable), 1);
TSFTable.SourceWeightSmoothness = repmat(iterative_smoothness, height(TSFTable), 1);
TSFTable.NoImprovePatience = repmat(iterative_no_improve_patience, height(TSFTable), 1);
TSFTable.TotalNLLFixedHypsometry = repmat(nll_fixed, height(TSFTable), 1);
TSFTable.TotalNLLSelectedTSF = repmat(nll_total, height(TSFTable), 1);
TSFTable.NLLChangeFromFixed = repmat(nll_total-nll_fixed, height(TSFTable), 1);
writetable(TSFTable, tsf_file);
fprintf("Wrote %s\n\n", tsf_file);

outFigDir = fullfile(catchment_dir, "figures_svg");
if ~exist(outFigDir, "dir"); mkdir(outFigDir); end
fig_tsf = figure('Name', 'Source_Weights');
subplot(1,2,1); hold on; box on;
if tsf_mode == "iterative"
    for g = 1:numel(tsf_group_names)
        member0 = tsf_group_indices{g}(1);
        plot(tsf_weights(:,member0) ./ max(pz, realmin), z_centers, ...
            'o-', 'LineWidth', 1.5, 'DisplayName', tsf_group_names(g));
    end
else
    plot(tsf_weights ./ max(pz, realmin), z_centers, ...
        'o-', 'LineWidth', 1.5, 'DisplayName', 'shared');
end
xline(1, '--k', 'DisplayName', 'hypsometry');
xlabel('Relative source weight (model / hypsometry)');
ylabel('Elevation (m)');
title(sprintf('Source weights (%s)', tsf_mode));
legend('Location', 'best');
subplot(1,2,2); hold on; box on;
stairs(cumsum(pz), z_edges(2:end), 'LineWidth', 1.5, ...
    'DisplayName', 'hypsometry');
if tsf_mode == "iterative"
    for g = 1:numel(tsf_group_names)
        member0 = tsf_group_indices{g}(1);
        stairs(cumsum(tsf_weights(:,member0)), z_edges(2:end), ...
            'LineWidth', 1.5, 'DisplayName', tsf_group_names(g));
    end
else
    stairs(cumsum(tsf_weights), z_edges(2:end), 'LineWidth', 1.5, ...
        'DisplayName', 'selected model TSF');
end
xlabel('Cumulative probability'); ylabel('Elevation (m)');
legend('Location', 'best');
title('Cumulative source weighting');
tsf_fig_file = fullfile(outFigDir, "Source_Weights.svg");
try
    exportgraphics(fig_tsf, tsf_fig_file, 'ContentType', 'vector');
catch
    print(fig_tsf, tsf_fig_file, '-dsvg');
end


%% ---------------- POST-FIT ORDERING DIAGNOSTICS ----------------
fprintf("Post-fit ordering check:\n");
PenaltyTable = table();
total_penalty = 0;
for p = 1:n_pairs
    lo = pair_lo(p);  hi = pair_hi(p);
    if ~(ok_mask(lo) && ok_mask(hi)); continue; end
    diff_vec  = Ahat_all(:,lo) - Ahat_all(:,hi);
    min_sep   = delta_min_Ma * pair_gap(p);
    viol      = max(0, diff_vec + min_sep);
    pen       = w_order * sum(viol.^2);
    total_penalty = total_penalty + pen;
    fprintf("  %s < %s  |  max violation=%.4f Ma  |  penalty=%.4f\n", ...
        chron_data(lo).chron, chron_data(hi).chron, max(viol), pen);
    row = table(string(chron_data(lo).chron), string(chron_data(hi).chron), ...
        pair_gap(p), delta_min_Ma*pair_gap(p), max(viol), pen, ...
        'VariableNames', {'Chron_Lo','Chron_Hi','Tc_Gap_Norm', ...
                          'MinSep_Ma','MaxViolation_Ma','PenaltyContribution'});
    PenaltyTable = [PenaltyTable; row]; %#ok<AGROW>
end
fprintf("  Total ordering penalty = %.4f\n\n", total_penalty);
pen_file = fullfile(catchment_dir, "ordering_penalty_contributions.csv");
writetable(PenaltyTable, pen_file);
fprintf("Wrote %s\n\n", pen_file);


%% ---------------- BOOTSTRAP ----------------
% Each resample draws new grains for ALL chronometers, then runs the full
% joint fminunc. This correctly propagates joint uncertainty.

Ab_all       = nan(Nz, n_chron, n_boot);
taub_all     = nan(n_chron, n_boot);
tsf_boot_all = nan(Nz, n_boot);
source_weights_boot = nan(Nz, n_chron, n_boot);
boot_tsf_fixed_nll    = nan(n_boot, 1);
boot_tsf_selected_nll = nan(n_boot, 1);
boot_tsf_best_iter    = nan(n_boot, 1);
boot_tsf_eval_iters   = nan(n_boot, 1);
boot_tsf_termination  = strings(n_boot, 1);
boot_error_message    = strings(n_boot, 1);
parallel_bootstrap_active = false;
bootstrap_workers_used    = 0;

if do_bootstrap
    fprintf("Phase 3: joint bootstrap (%d resamples)...\n", n_boot);
    fprintf("  Random seed: %d\n", bootstrap_random_seed);
    fprintf("  (Each resample is a full joint fminunc solve -- may take a few minutes)\n");
    if tsf_mode == "iterative"
        fprintf("  Configured source-weight groups are re-estimated in every resample.\n");
        fprintf("  Iterative weight updates are capped at %d per resample.\n", ...
            iterative_boot_max_outer);
    end

    opts_boot = optimoptions('fminunc', ...
        'Algorithm',              'quasi-newton', ...
        'MaxIterations',          3000, ...
        'MaxFunctionEvaluations', 5e5, ...
        'FiniteDifferenceType',   'central', ...
        'Display',                'off', ...
        'OptimalityTolerance',    1e-7, ...
        'StepTolerance',          1e-9);

    if use_parallel_bootstrap
        has_parallel_tools = exist('parpool', 'file') == 2 && ...
            license('test', 'Distrib_Computing_Toolbox');
        if has_parallel_tools
            try
                can_start_pool = true;
                if exist('canUseParallelPool', 'file') == 2
                    can_start_pool = canUseParallelPool;
                end
                if ~can_start_pool
                    error("A local parallel pool is not currently available.");
                end
                pool = gcp('nocreate');
                if isempty(pool)
                    if isempty(parallel_worker_count)
                        pool = parpool('Processes');
                    else
                        pool = parpool('Processes', parallel_worker_count);
                    end
                elseif ~isempty(parallel_worker_count) && ...
                        pool.NumWorkers ~= parallel_worker_count
                    warning("An existing pool has %d workers; using it instead of the requested %d.", ...
                        pool.NumWorkers, parallel_worker_count);
                end
                parallel_bootstrap_active = pool.NumWorkers > 1;
                bootstrap_workers_used = pool.NumWorkers;
            catch ME
                warning("Parallel bootstrap unavailable (%s). Continuing serially.", ...
                    ME.message);
            end
        else
            warning("Parallel Computing Toolbox is unavailable. Continuing serially.");
        end
    end

    report_bootstrap_progress(n_boot, true);
    if parallel_bootstrap_active
        fprintf("  Parallel execution: %d process workers\n", ...
            bootstrap_workers_used);
        progress_queue = parallel.pool.DataQueue;
        progress_listener = afterEach(progress_queue, ...
            @(~) report_bootstrap_progress(n_boot, false));
        parfor b = 1:n_boot
            [Ab_b, taub_b, tsf_vector_b, tsf_matrix_b, fixed_nll_b, ...
                selected_nll_b, best_iter_b, eval_iters_b, termination_b, ...
                error_b] = run_bootstrap_replicate( ...
                    b, bootstrap_random_seed, chron_data, ok_idx, pz, Nz, ...
                    pair_lo, pair_hi, pair_gap, w_order, delta_min_Ma, ...
                    theta_fixed_hat, opts_boot, tsf_mode, ...
                    tsf_group_indices, tsf_group_fixed, age_grid_step, ...
                    iterative_hypsometry_pull, iterative_smoothness, ...
                    tsf_update_fraction, tsf_smooth_span, ...
                    iterative_boot_max_outer, iterative_weight_tol, ...
                    iterative_nll_tol, iterative_no_improve_patience);
            Ab_all(:,:,b)              = Ab_b;
            taub_all(:,b)              = taub_b;
            tsf_boot_all(:,b)          = tsf_vector_b;
            source_weights_boot(:,:,b) = tsf_matrix_b;
            boot_tsf_fixed_nll(b)      = fixed_nll_b;
            boot_tsf_selected_nll(b)   = selected_nll_b;
            boot_tsf_best_iter(b)      = best_iter_b;
            boot_tsf_eval_iters(b)     = eval_iters_b;
            boot_tsf_termination(b)    = termination_b;
            boot_error_message(b)      = error_b;
            send(progress_queue, b);
        end
    else
        fprintf("  Serial execution\n");
        for b = 1:n_boot
            [Ab_all(:,:,b), taub_all(:,b), tsf_boot_all(:,b), ...
                source_weights_boot(:,:,b), boot_tsf_fixed_nll(b), ...
                boot_tsf_selected_nll(b), boot_tsf_best_iter(b), ...
                boot_tsf_eval_iters(b), boot_tsf_termination(b), ...
                boot_error_message(b)] = run_bootstrap_replicate( ...
                    b, bootstrap_random_seed, chron_data, ok_idx, pz, Nz, ...
                    pair_lo, pair_hi, pair_gap, w_order, delta_min_Ma, ...
                    theta_fixed_hat, opts_boot, tsf_mode, ...
                    tsf_group_indices, tsf_group_fixed, age_grid_step, ...
                    iterative_hypsometry_pull, iterative_smoothness, ...
                    tsf_update_fraction, tsf_smooth_span, ...
                    iterative_boot_max_outer, iterative_weight_tol, ...
                    iterative_nll_tol, iterative_no_improve_patience);
            report_bootstrap_progress(n_boot, false);
        end
    end

    for b = 1:n_boot
        if strlength(boot_error_message(b)) > 0
            fprintf("  Bootstrap %d failed: %s\n", b, boot_error_message(b));
        end
    end
    valid_boot = squeeze(~any(any(isnan(Ab_all(:,ok_idx,:)), 1), 2));
    n_boot_success = sum(valid_boot);
    fprintf("  Bootstrap complete: %d / %d successful resamples.\n\n", ...
        n_boot_success, n_boot);

    if tsf_mode == "iterative"
        boot_success = valid_boot(:);
        BootTSFConvergence = table( ...
            (1:n_boot)', boot_success, boot_tsf_fixed_nll, ...
            boot_tsf_selected_nll, ...
            boot_tsf_selected_nll - boot_tsf_fixed_nll, ...
            boot_tsf_best_iter, boot_tsf_eval_iters, boot_tsf_termination, ...
            'VariableNames', {'BootstrapReplicate','Successful', ...
                'FixedHypsometryNLL','SelectedIterativeNLL', ...
                'NLLChangeFromFixed','BestIteration', ...
                'EvaluatedIterations','TerminationReason'});
        boot_conv_file = fullfile( ...
            catchment_dir, "source_weighting_bootstrap_summary.csv");
        writetable(BootTSFConvergence, boot_conv_file);
        fprintf("Wrote %s\n\n", boot_conv_file);
    end
end

if do_bootstrap && tsf_mode == "iterative"
    SourceWeightBootTable = table();
    for g = 1:numel(tsf_group_names)
        member0 = tsf_group_indices{g}(1);
        boot_g = reshape(source_weights_boot(:,member0,:), Nz, n_boot);
        valid_g = boot_g(:, ~any(isnan(boot_g), 1));
        if isempty(valid_g)
            continue;
        end

        boot_median_g = median(valid_g, 2, 'omitnan');
        boot_lo_g     = quantile(valid_g, ci_lo, 2);
        boot_hi_g     = quantile(valid_g, ci_hi, 2);
        group_name_col = repmat(tsf_group_names(g), Nz, 1);
        fixed_col = repmat(tsf_group_fixed(g), Nz, 1);
        n_success_col = repmat(size(valid_g, 2), Nz, 1);

        GroupWeightBoot = table( ...
            group_name_col, fixed_col, z_centers, pz, ...
            tsf_weights(:,member0), boot_median_g, boot_lo_g, boot_hi_g, ...
            n_success_col, ...
            'VariableNames', {'TSFGroup','FixedToHypsometry', ...
                'Elevation_m','HypsometryWeight','BestFitTSFWeight', ...
                'TSF_boot_median','TSF_boot_p16','TSF_boot_p84', ...
                'SuccessfulResamples'});
        SourceWeightBootTable = [SourceWeightBootTable; GroupWeightBoot]; %#ok<AGROW>
    end

    if ~isempty(SourceWeightBootTable)
        source_weight_boot_file = fullfile( ...
            catchment_dir, "source_weighting_bootstrap_CI.csv");
        writetable(SourceWeightBootTable, source_weight_boot_file);
        fprintf("Wrote %s\n\n", source_weight_boot_file);
    end
end
%% ---------------- PER-CHRONOMETER OUTPUTS ----------------
summary_file = fullfile(catchment_dir, "summary_fit_params.csv");
Summary      = table();

opts_disp = optimset('MaxIter', 1e4, 'MaxFunEvals', 1e5, 'Display', 'off');

for c = ok_idx

    cd     = chron_data(c);
    chron  = cd.chron;
    age    = cd.age;
    sig    = cd.sig;
    N      = cd.N;
    Ahat   = Ahat_all(:, c);
    tauhat = tauhat_all(c);
    tsf_weights_c = tsf_for_chron(tsf_weights, c);

    fprintf("===== %s  (joint fit) =====\n", chron);
    fprintf("  tau = %.3f Ma  |  A = [%.2f, %.2f] Ma\n", ...
        tauhat, min(Ahat), max(Ahat));

    % Per-chronometer likelihood contributions for selected and reference fits
    nll_c = nll_single(theta_joint_hat( block_slice(c, ok_idx, Nz) ), ...
                       age, sig, tsf_weights_c, cd.AminB, cd.AmaxB, ...
                       cd.tau_min, cd.lambda);
    nll_c_fixed = nll_single(theta_fixed_hat( block_slice(c, ok_idx, Nz) ), ...
                       age, sig, pz, cd.AminB, cd.AmaxB, ...
                       cd.tau_min, cd.lambda);
    fprintf("  Single-chron NLL contribution = %.4f", nll_c);
    if tsf_mode ~= "fixed"
        fprintf("  (fixed reference %.4f)", nll_c_fixed);
    end
    fprintf("\n\n");

    % Track figures created in this block
    figs_before = findall(0, 'Type', 'figure');
    nums_before = arrayfun(@(h) h.Number, figs_before);

    %% ---- Bootstrap CI ----
    if do_bootstrap
        Ab_c  = squeeze(Ab_all(:, c, :));   % Nz x n_boot
        Ab_c_valid = Ab_c(:, ~any(isnan(Ab_c), 1));

        A_med = median(Ab_c_valid, 2, 'omitnan');
        A_lo  = quantile(Ab_c_valid, ci_lo, 2);
        A_hi  = quantile(Ab_c_valid, ci_hi, 2);

        taub_c   = squeeze(taub_all(c, :));
        tau_med  = median(taub_c, 'omitnan');
        tau_lo   = quantile(taub_c, ci_lo);
        tau_hi   = quantile(taub_c, ci_hi);

        fprintf("  Bootstrap tau (median [%.0f%% CI]): %.3f [%.3f, %.3f]\n", ...
            (ci_hi-ci_lo)*100, tau_med, tau_lo, tau_hi);

        A_sigma_ish   = 0.5 * (A_hi - A_lo);
        A_sigma_minus = A_med - A_lo;
        A_sigma_plus  = A_hi  - A_med;

        out_ci = table(z_centers, pz, tsf_weights_c, Ahat_fixed(:,c), Ahat, ...
            A_med, A_lo, A_hi, ...
            A_sigma_ish, A_sigma_minus, A_sigma_plus, ...
            'VariableNames', {'Elevation_m','HypsometryWeight','ModelTSFWeight', ...
                              'Ahat_fixed_hypsometry_Ma','Ahat_Ma', ...
                              'A_boot_median_Ma','A_boot_p16_Ma','A_boot_p84_Ma', ...
                              'A_sigma_ish_Ma','A_sigma_minus_Ma','A_sigma_plus_Ma'});
        out_ci.A_med_pm1sig = arrayfun( ...
            @(m,s) sprintf('%.2f +/- %.2f', m, s), A_med, A_sigma_ish, ...
            'UniformOutput', false);
        out_ci.A_med_asym1sig = arrayfun( ...
            @(m,sp,sm) sprintf('%.2f +%.2f/-%.2f', m, sp, sm), ...
            A_med, A_sigma_plus, A_sigma_minus, 'UniformOutput', false);

        ci_file = fullfile(catchment_dir, ...
            "predicted_bedrock_transect_" + string(chron) + "_CI.csv");
        writetable(out_ci, ci_file);
        fprintf("  Wrote %s\n", ci_file);

        % Bootstrap plot
        figure('Name', sprintf('%s_Az_bootstrap', chron)); hold on;
        fill([z_centers; flipud(z_centers)], [A_lo; flipud(A_hi)], ...
            [0.8 0.8 0.8], 'EdgeColor', 'none');
        plot(z_centers, Ahat, 'k-o', 'LineWidth', 1.5);
        xlabel('Elevation (m)'); ylabel('Predicted bedrock age (Ma)');
        title(sprintf('%s: A(z) with joint bootstrap CI', chron));
        legend('68% CI', 'Best-fit A(z)', 'Location', 'best');
    end

    %% ---- Predicted vs observed detrital PDF ----
    amin     = max(0, floor(min(age) - 5));
    amax     = ceil(max(age) + 5);
    age_grid = (amin : age_grid_step : amax)';

    y_obs = zeros(size(age_grid));
    for i = 1:N
        y_obs = y_obs + normpdf(age_grid, age(i), max(sig(i), 1e-6));
    end
    y_obs = y_obs / trapz(age_grid, y_obs);

    sig_pred = sqrt(median(sig)^2 + tauhat^2);
    y_pred   = zeros(size(age_grid));
    for k = 1:Nz
        y_pred = y_pred + tsf_weights_c(k) * ...
            normpdf(age_grid, Ahat(k), sig_pred);
    end
    y_pred = y_pred / trapz(age_grid, y_pred);

    figure('Name', sprintf('%s_PDF_Obs_vs_Pred', chron)); hold on;
    plot(age_grid, y_obs, 'LineWidth', 1.5);
    plot(age_grid, y_pred, 'LineWidth', 1.5);
    xlabel('Age (Ma)'); ylabel('Probability density');
    legend('Observed (kernel)', 'Predicted (selected TSF)', 'Location', 'best');
    title(sprintf('%s: observed vs predicted detrital PDF', chron));

    %% ---- Per-grain posterior p(z_k | age_i) ----
    Pzik = zeros(Nz, N);
    for i = 1:N
        s      = sqrt(sig(i)^2 + tauhat^2);
        like_k = normpdf(age(i), Ahat, s);
        post   = like_k(:) .* tsf_weights_c(:);
        post   = post / max(sum(post), realmin);
        Pzik(:, i) = post;
    end
    z_mean_per_grain = (z_centers(:)' * Pzik)';
    pz_implied       = mean(Pzik, 2);
    pz_implied       = pz_implied / sum(pz_implied);

    z_p05 = zeros(N,1); z_p16 = zeros(N,1); z_p50 = zeros(N,1);
    z_p84 = zeros(N,1); z_p95 = zeros(N,1); z_sd  = zeros(N,1);
    for i = 1:N
        p   = Pzik(:,i);
        cdf = cumsum(p);
        [cdfu, ia_i] = unique(cdf, 'stable');
        zu_i = z_centers(ia_i);
        if numel(cdfu) < 2
            zm = sum(p .* z_centers);
            z_p05(i)=zm; z_p16(i)=zm; z_p50(i)=zm;
            z_p84(i)=zm; z_p95(i)=zm; z_sd(i)=0;
            continue;
        end
        z_p05(i) = interp1(cdfu, zu_i, 0.05, 'linear', 'extrap');
        z_p16(i) = interp1(cdfu, zu_i, 0.16, 'linear', 'extrap');
        z_p50(i) = interp1(cdfu, zu_i, 0.50, 'linear', 'extrap');
        z_p84(i) = interp1(cdfu, zu_i, 0.84, 'linear', 'extrap');
        z_p95(i) = interp1(cdfu, zu_i, 0.95, 'linear', 'extrap');
        zm = sum(p .* z_centers);
        z_sd(i) = sqrt(sum(p .* (z_centers - zm).^2));
    end

    %% ---- Remaining diagnostic figures ----
    figure('Name', sprintf('%s_Az_curve', chron)); hold on;
    plot(z_centers, Ahat, 'k-o', 'LineWidth', 1.5);
    xlabel('Elevation (m)'); ylabel('Predicted bedrock age A(z_k) (Ma)');
    title(sprintf('%s: joint-fit age-elevation curve (%s)', chron, tsf_mode));

    figure('Name', sprintf('%s_CDF_Hyps_vs_Implied', chron)); hold on;
    stairs(z_edges(2:end), cumsum(pz),         'LineWidth', 1.5);
    if tsf_mode ~= "fixed"
        stairs(z_edges(2:end), cumsum(tsf_weights_c), 'LineWidth', 1.5);
    end
    stairs(z_edges(2:end), cumsum(pz_implied), 'LineWidth', 1.5);
    xlabel('Elevation (m)'); ylabel('Cumulative probability');
    if tsf_mode ~= "fixed"
        legend('Hypsometry CDF', 'Selected model TSF', ...
               'Post-fit posterior diagnostic', 'Location', 'best');
    else
        legend('Hypsometry CDF', 'Posterior diagnostic', 'Location', 'best');
    end
    title(sprintf('%s: source-weighting diagnostics', chron));

    figure('Name', sprintf('%s_Grain_Ez', chron));
    scatter(age, z_mean_per_grain, 20, 'filled');
    xlabel('Detrital age (Ma)'); ylabel('E[source elevation | age] (m)');
    title(sprintf('%s: grain-wise expected source elevation', chron));

    %% ---- Save figures ----
    figs_after = findall(0, 'Type', 'figure');
    nums_after = arrayfun(@(h) h.Number, figs_after);
    new_nums   = setdiff(nums_after, nums_before);

    outFigDir = fullfile(catchment_dir, "figures_svg");
    if ~exist(outFigDir, "dir"); mkdir(outFigDir); end

    for nn = 1:numel(new_nums)
        fig     = figure(new_nums(nn));
        figname = string(get(fig, 'Name'));
        if strlength(figname) == 0; figname = "Figure"; end
        figname = regexprep(figname, '[^\w\-]+', '_');
        outfile = fullfile(outFigDir, ...
            sprintf("%s_fig%02d_%s.svg", chron, fig.Number, figname));
        try
            exportgraphics(fig, outfile, 'ContentType', 'vector');
        catch
            print(fig, outfile, '-dsvg');
        end
    end
    fprintf("  Saved %d figures to %s\n", numel(new_nums), outFigDir);

    %% ---- Write CSV outputs ----
    out_transect = table(z_centers, pz, tsf_weights_c, Ahat, Ahat_fixed(:,c), ...
        'VariableNames', {'Elevation_m','HypsometryWeight','ModelTSFWeight', ...
                          'PredictedBedrockAge_Ma', ...
                          'FixedHypsometryBedrockAge_Ma'});
    writetable(out_transect, fullfile(catchment_dir, ...
        "predicted_bedrock_transect_" + string(chron) + ".csv"));
    fprintf("  Wrote transect CSV\n");

    if write_grain_posteriors
        Pfile  = fullfile(catchment_dir, "grain_posteriors_" + string(chron) + ".csv");
        Ptable = array2table(Pzik);
        Ptable = addvars(Ptable, z_centers, 'Before', 1, 'NewVariableNames', "Elevation_m");
        Ptable = addvars(Ptable, tsf_weights_c, 'After', "Elevation_m", ...
            'NewVariableNames', "ModelTSFWeight");
        writetable(Ptable, Pfile);

        GrainT = table(age, sig, z_mean_per_grain, z_p50, z_p16, z_p84, z_p05, z_p95, z_sd, ...
            'VariableNames', {'Date_Ma','Error_Ma','E_SourceElev_m','Median_SourceElev_m', ...
                              'P16_SourceElev_m','P84_SourceElev_m','P05_SourceElev_m', ...
                              'P95_SourceElev_m','SD_SourceElev_m'});
        writetable(GrainT, fullfile(catchment_dir, ...
            "grain_expected_source_" + string(chron) + ".csv"));
        fprintf("  Wrote grain posterior CSVs\n");
    end

    % Append to summary
    chron_str = string(chron);
    tsf_tv = 0.5 * sum(abs(tsf_weights_c - pz));
    Summary = [Summary; table(chron_str, N, cd.n_excluded, ...
        nll_c, nll_c_fixed, nll_c-nll_c_fixed, tauhat, tauhat_fixed(c), ...
        min(Ahat), max(Ahat), tsf_mode, chron_tsf_group(c), tsf_tv, ...
        cd.tau_min, cd.age_margin, cd.lambda, ...
        cd.age_min_filt, cd.age_max_filt, w_order, delta_min_Ma, Tc_vec(c), ...
        tsf_update_fraction, tsf_smooth_span, ...
        iterative_hypsometry_pull, iterative_smoothness, ...
        iterative_no_improve_patience, do_bootstrap, n_boot, ...
        bootstrap_random_seed, use_parallel_bootstrap, ...
        bootstrap_workers_used, ...
        'VariableNames', {'Chronometer','Ngrains','Nexcluded', ...
                          'NLL_single','NLL_single_fixed_hypsometry', ...
                          'NLL_change_from_fixed','Tau_Ma', ...
                          'Tau_fixed_hypsometry_Ma','Amin_Ma','Amax_Ma', ...
                          'TSFMode','TSFGroup', ...
                          'TSF_TotalVariationFromHypsometry', ...
                          'Setting_TauMin', ...
                          'Setting_AgeMargin','Setting_Lambda', ...
                          'Setting_AgeMinFilter','Setting_AgeMaxFilter', ...
                          'Setting_w_order','Setting_delta_min_Ma','Tc_degC', ...
                          'Setting_SourceWeightUpdateFraction', ...
                          'Setting_SourceWeightSmoothSpan', ...
                          'Setting_IterativeHypsometryPull', ...
                          'Setting_IterativeSmoothness', ...
                          'Setting_IterativeNoImprovePatience', ...
                          'Setting_DoBootstrap', 'Setting_NBootstrap', ...
                          'Setting_BootstrapRandomSeed', ...
                          'Setting_ParallelBootstrapRequested', ...
                          'Setting_ParallelWorkersUsed'})]; %#ok<AGROW>

    fprintf("\n");
end


%% ---------------- JOINT ORDERING SUMMARY FIGURE ----------------
% Plot all chronometer A(z) curves together to visualise ordering.
if n_ok >= 2
    colors_map = {'#378ADD','#E24B4A','#1D9E75','#BA7517', ...
                  '#7F77DD','#D85A30','#0F6E56','#633806'};
    figure('Name', 'AllChrons_AgeElevation_Joint'); hold on; box on; grid on;
    leg_entries = {};
    ci_idx = 0;
    for c = ok_idx
        ci_idx = ci_idx + 1;
        col = colors_map{mod(ci_idx-1, numel(colors_map))+1};
        plot(Ahat_all(:,c), z_centers, 'o-', 'Color', col, ...
            'LineWidth', 1.5, 'MarkerFaceColor', col, 'MarkerSize', 5);
        leg_entries{end+1} = chron_data(c).chron; %#ok<AGROW>
    end
    xlabel('Age (Ma)'); ylabel('Elevation (m)');
    title(sprintf('%s: all chronometers — %s source weighting', ...
        catchment_name, tsf_mode));
    legend(leg_entries, 'Location', 'best');

    outFigDir = fullfile(catchment_dir, "figures_svg");
    if ~exist(outFigDir, "dir"); mkdir(outFigDir); end
    outfile = fullfile(outFigDir, "AllChrons_AgeElevation_Joint.svg");
    try
        exportgraphics(gcf, outfile, 'ContentType', 'vector');
    catch
        print(gcf, outfile, '-dsvg');
    end
end


%% ---------------- SAVE SUMMARY ----------------
if ~isempty(Summary)
    writetable(Summary, summary_file);
    fprintf("Wrote summary : %s\n", summary_file);
else
    fprintf("No chronometers processed.\n");
end

fprintf("\nDone.\n");


%% ===============================================================
%% LOCAL FUNCTIONS
%% ===============================================================

function nll = nll_single(theta, age, sig, pz, AminB, AmaxB, tau_min, lambda)
% NLL for a single chronometer. Grain-by-bin likelihoods are evaluated as
% one matrix because this function is called many times by the optimizer.
    Nz = numel(pz);
    [A, tau] = unpack_single(theta, Nz, AminB, AmaxB, tau_min);
    age = age(:)';
    combined_sigma = sqrt(sig(:)'.^2 + tau^2);
    density = gaussian_density(age, A(:), combined_sigma);
    mixture_probability = sum(pz(:) .* density, 1);
    nll = -sum(log(max(mixture_probability, realmin)));
    if lambda > 0 && Nz >= 3
        d2  = A(3:end) - 2*A(2:end-1) + A(1:end-2);
        nll = nll + lambda * sum(d2.^2);
    end
end


function [A, tau] = unpack_single(theta, Nz, AminB, AmaxB, tau_min)
% Decode single-chronometer parameter vector into bounded monotonic A(z).
    s          = theta(1:Nz);
    logtau_raw = theta(Nz+1);
    Araw = AminB + (AmaxB - AminB) .* (1 ./ (1 + exp(-s)));
    A    = Araw;
    for k = 2:Nz
        A(k) = max(A(k), A(k-1));
    end
    A(A > AmaxB) = AmaxB;
    A(A < AminB) = AminB;
    tau = tau_min + exp(logtau_raw);
end


function nll = nll_joint(theta, chron_data, ok_idx, tsf_weights, Nz, ...
                          pair_lo, pair_hi, pair_gap, w_order, delta_min_Ma)
% Joint NLL: sum of per-chronometer NLLs + pairwise ordering penalty.
% tsf_weights may be one shared Nz-vector or an Nz-by-n_chron matrix.

    % Unpack all A(z) curves
    A_all   = nan(Nz, numel(chron_data));
    nll     = 0;
    pos     = 1;

    for c = ok_idx
        cd    = chron_data(c);
        blk   = pos : pos + Nz;           % Nz+1 elements
        th_c  = theta(blk);
        pos   = pos + Nz + 1;

        [A_c, ~] = unpack_single(th_c, Nz, cd.AminB, cd.AmaxB, cd.tau_min);
        A_all(:, c) = A_c;

        w_c = tsf_for_chron(tsf_weights, c);
        nll = nll + nll_single(th_c, cd.age, cd.sig, w_c, ...
                               cd.AminB, cd.AmaxB, cd.tau_min, cd.lambda);
    end

    % Pairwise ordering penalty
    for p = 1:numel(pair_lo)
        lo = pair_lo(p);  hi = pair_hi(p);
        if any(isnan(A_all(:,lo))) || any(isnan(A_all(:,hi))); continue; end
        min_sep  = delta_min_Ma * pair_gap(p);
        viol     = max(0, A_all(:,lo) - A_all(:,hi) + min_sep);
        nll      = nll + w_order * sum(viol.^2);
    end
end


function w_c = tsf_for_chron(tsf_weights, c)
% Return the source-weight vector used by chronometer c.
    if isvector(tsf_weights)
        w_c = tsf_weights(:);
    else
        w_c = tsf_weights(:,c);
    end
    w_c = w_c(:);
end


function [group_names, group_indices, group_fixed, chron_group] = ...
        build_configured_tsf_groups(chron_data, ok_idx)
% Build effective source-weight groups from explicit config values. This
% avoids assigning scientific meaning from a chronometer's display name.

    n_chron = numel(chron_data);
    chron_group = repmat("unassigned", n_chron, 1);
    labels = strings(numel(ok_idx), 1);
    keys = strings(numel(ok_idx), 1);
    estimate = false(numel(ok_idx), 1);

    for j = 1:numel(ok_idx)
        c = ok_idx(j);
        labels(j) = strtrim(string(chron_data(c).tsf_group));
        keys(j) = lower(labels(j));
        estimate(j) = chron_data(c).estimate_tsf;
    end

    [unique_keys, first_idx] = unique(keys, 'stable');
    group_names = labels(first_idx);
    group_indices = cell(numel(unique_keys), 1);
    group_fixed = false(numel(unique_keys), 1);

    for g = 1:numel(unique_keys)
        in_group = keys == unique_keys(g);
        member_positions = find(in_group);
        members = ok_idx(member_positions);
        flags = estimate(in_group);
        if any(flags ~= flags(1))
            error("All chronometers in TSFGroup '%s' must use the same EstimateTSF value.", ...
                group_names(g));
        end
        group_indices{g} = members;
        group_fixed(g) = ~flags(1);
        chron_group(members) = group_names(g);
    end
end


function value = parse_config_logical(raw_value, column_name, chronometer)
% Parse true/false config values without depending on readtable's inferred
% column type. Accepted forms: true/false, yes/no, 1/0.

    if islogical(raw_value)
        value = raw_value;
        return;
    end
    if isnumeric(raw_value)
        if isscalar(raw_value) && isfinite(raw_value) && any(raw_value == [0, 1])
            value = logical(raw_value);
            return;
        end
    end

    text_value = lower(strtrim(string(raw_value)));
    if ismissing(text_value) || strlength(text_value) == 0
        error("%s is required for chronometer '%s' in iterative mode.", ...
            column_name, chronometer);
    elseif any(text_value == ["true", "yes", "1"])
        value = true;
    elseif any(text_value == ["false", "no", "0"])
        value = false;
    else
        error("%s for chronometer '%s' must be true/false, yes/no, or 1/0.", ...
            column_name, chronometer);
    end
end


function [A_all, tau_all] = unpack_joint(theta, chron_data, ok_idx, Nz)
% Unpack full joint theta into per-chronometer A(z) and tau arrays.
    n_chron = numel(chron_data);
    A_all   = nan(Nz, n_chron);
    tau_all = nan(n_chron, 1);
    pos     = 1;
    for c = ok_idx
        cd  = chron_data(c);
        blk = pos : pos + Nz;
        [A_all(:,c), tau_all(c)] = unpack_single(theta(blk), Nz, ...
                                       cd.AminB, cd.AmaxB, cd.tau_min);
        pos = pos + Nz + 1;
    end
end


function idx = block_slice(c, ok_idx, Nz)
% Return the index range in the joint theta vector for chronometer c.
    pos = 1;
    for ci = ok_idx
        if ci == c
            idx = pos : pos + Nz;
            return;
        end
        pos = pos + Nz + 1;
    end
    idx = [];
end


function [w_update, n_grains] = posterior_tsf_update( ...
        A_all, tau_all, chron_data, ok_idx, w_prior)
% One pooled EM-style update of shared elevation-mixture weights.
%
% Responsibilities are evaluated under the supplied A(z), tau, and prior
% weights, then averaged over all retained grains and chronometers. Pooling
% by grain matches the weighting used by the joint log likelihood, so a
% chronometer with more retained grains contributes proportionally more.

    Nz          = numel(w_prior);
    weight_sum  = zeros(Nz, 1);
    n_grains    = 0;

    for c = ok_idx
        cd = chron_data(c);
        combined_sigma = sqrt(cd.sig(:)'.^2 + tau_all(c)^2);
        likelihood = gaussian_density( ...
            cd.age(:)', A_all(:,c), combined_sigma);
        weighted_likelihood = w_prior(:) .* likelihood;
        responsibility = weighted_likelihood ./ ...
            max(sum(weighted_likelihood, 1), realmin);
        weight_sum = weight_sum + sum(responsibility, 2);
        n_grains   = n_grains + cd.N;
    end

    if n_grains == 0 || sum(weight_sum) <= 0
        w_update = w_prior(:) / sum(w_prior);
        return;
    end

    w_update = weight_sum / n_grains;
    w_update = max(w_update, 0);
    w_update = w_update / sum(w_update);
end


function [estimated_weights, data_sse, n_data_rows] = estimate_source_weights_nnls( ...
        A_all, tau_all, chron_data, ok_idx, w_hypso, age_step, ...
        hypsometry_pull, smoothness)
% Estimate one source-weight curve with regularized nonnegative least squares.
%
% For each chronometer, columns of G are predicted age densities from the
% individual elevation bins and y is the observed detrital age density.
% The blocks are scaled by grain count and normalized density magnitude,
% then stacked. The augmented least-squares system applies:
%   - nonnegativity through lsqnonneg;
%   - a quadratic pull toward hypsometry;
%   - a second-difference smoothness penalty;
%   - an approximate sum-to-one constraint, followed by exact normalization.

    Nz = numel(w_hypso);
    total_grains = sum(arrayfun(@(c) chron_data(c).N, ok_idx));
    G_all = [];
    y_all = [];

    for c = ok_idx
        cd  = chron_data(c);
        age = cd.age(:);
        sig = cd.sig(:);
        A   = A_all(:,c);
        tau = tau_all(c);

        max_sd = max(sqrt(sig.^2 + tau.^2));
        amin = max(0, min([age; A]) - 4 * max_sd);
        amax = max([age; A]) + 4 * max_sd;
        if amax <= amin
            amax = amin + max(age_step, 1);
        end
        grid = (amin:age_step:amax)';
        if grid(end) < amax
            grid(end+1,1) = amax;
        end

        % Observed KDE retains each grain's reported analytical error.
        y = sum(gaussian_density( ...
            grid, age', max(sig', 1e-6)), 2);
        y_area = trapz(grid, y);
        if y_area <= 0 || ~isfinite(y_area)
            continue;
        end
        y = y / y_area;

        % Predicted density for each elevation averages over the analytical
        % error distribution of the grains and includes fitted extra scatter.
        G = zeros(numel(grid), Nz);
        combined_sigma = sqrt(sig'.^2 + tau^2);
        for k = 1:Nz
            G(:,k) = mean(gaussian_density( ...
                grid, A(k), combined_sigma), 2);
            g_area = trapz(grid, G(:,k));
            if g_area > 0 && isfinite(g_area)
                G(:,k) = G(:,k) / g_area;
            end
        end

        % Equalize density units while retaining grain-count weighting.
        density_norm = max(norm(sqrt(age_step) * y), realmin);
        block_scale = sqrt(cd.N / max(total_grains, 1)) / density_norm;
        G_all = [G_all; block_scale * sqrt(age_step) * G]; %#ok<AGROW>
        y_all = [y_all; block_scale * sqrt(age_step) * y]; %#ok<AGROW>
    end

    if isempty(G_all)
        estimated_weights = w_hypso(:) / sum(w_hypso);
        data_sse = NaN;
        n_data_rows = 0;
        return;
    end

    data_scale = max(mean(sum(G_all.^2, 1)), realmin);
    C = G_all;
    d = y_all;

    if hypsometry_pull > 0
        pull_scale = sqrt(hypsometry_pull * data_scale);
        C = [C; pull_scale * eye(Nz)];
        d = [d; pull_scale * w_hypso(:)];
    end

    if smoothness > 0 && Nz >= 3
        D2 = diff(eye(Nz), 2, 1);
        smooth_scale = sqrt(smoothness * data_scale);
        C = [C; smooth_scale * D2];
        d = [d; zeros(Nz-2, 1)];
    end

    % Enforce sum(w)=1 strongly in the augmented system. Exact
    % normalization below removes the small residual allowed by this row.
    sum_scale = 100 * sqrt(data_scale);
    C = [C; sum_scale * ones(1, Nz)];
    d = [d; sum_scale];

    estimated_weights = lsqnonneg(C, d);
    if sum(estimated_weights) <= 0 || any(~isfinite(estimated_weights))
        estimated_weights = w_hypso(:);
    end
    estimated_weights = max(estimated_weights(:), 1e-12);
    estimated_weights = estimated_weights / sum(estimated_weights);

    data_sse = sum((G_all * estimated_weights - y_all).^2);
    n_data_rows = size(G_all, 1);
end


function density = gaussian_density(x, mu, sigma)
% Normal probability density with implicit expansion. Callers arrange x,
% mu, and sigma as row/column vectors to construct the required matrix.
    density = normpdf(x, mu, sigma);
end


function [theta_best, weights_best, best_nll, best_iteration, ...
        evaluated_iterations, termination_reason] = ...
        fit_iterative_source_weights_bootstrap( ...
            theta_fixed, nll_fixed, chron_data, ok_idx, pz, Nz, ...
            pair_lo, pair_hi, pair_gap, w_order, delta_min_Ma, opts_fit, ...
            group_indices, group_fixed, age_grid_step, hypsometry_pull, ...
            smoothness, update_fraction, smooth_span, max_outer, ...
            weight_tol, nll_tol, no_improve_patience)
% Guarded iterative source-weight inversion for one bootstrap resample.
% Each resample starts from its independently fitted fixed-hypsometry state.
% The source-weight density objective and grain likelihood are not identical, so the
% best-likelihood state is retained if later alternating updates deteriorate.

    n_chron = numel(chron_data);
    weights_current = repmat(pz(:), 1, n_chron);
    theta_current = theta_fixed;
    nll_previous = nll_fixed;

    theta_best    = theta_fixed;
    weights_best  = weights_current;
    best_nll      = nll_fixed;
    best_iteration = 0;
    no_improve_count = 0;
    converged = false;
    stopped_no_improve = false;
    evaluated_iterations = 0;

    for outer = 1:max_outer
        [A_current, tau_current] = unpack_joint( ...
            theta_current, chron_data, ok_idx, Nz);
        weights_new = weights_current;

        for g = 1:numel(group_indices)
            members = group_indices{g};
            if group_fixed(g)
                target_g = pz(:);
            else
                target_g = estimate_source_weights_nnls( ...
                    A_current, tau_current, chron_data, members, pz, ...
                    age_grid_step, hypsometry_pull, smoothness);
                if smooth_span > 1
                    target_g = smoothdata(target_g, 'movmean', smooth_span);
                end
                target_g = max(target_g, 1e-8);
                target_g = target_g / sum(target_g);
            end

            for c_group = members
                if group_fixed(g)
                    weights_new(:,c_group) = pz(:);
                else
                    weights_new(:,c_group) = ...
                        (1 - update_fraction) * weights_current(:,c_group) + ...
                        update_fraction * target_g;
                    weights_new(:,c_group) = max( ...
                        weights_new(:,c_group), 1e-8);
                    weights_new(:,c_group) = weights_new(:,c_group) / ...
                        sum(weights_new(:,c_group));
                end
            end
        end

        obj_iterative = @(th) nll_joint( ...
            th, chron_data, ok_idx, weights_new, Nz, ...
            pair_lo, pair_hi, pair_gap, w_order, delta_min_Ma);
        [theta_new, nll_new] = fminunc( ...
            obj_iterative, theta_current, opts_fit);

        max_dw = max(abs(weights_new(:) - weights_current(:)));
        nll_change = nll_new - nll_previous;
        evaluated_iterations = outer;

        if nll_new < best_nll
            theta_best       = theta_new;
            weights_best     = weights_new;
            best_nll         = nll_new;
            best_iteration   = outer;
            no_improve_count = 0;
        else
            no_improve_count = no_improve_count + 1;
        end

        theta_current   = theta_new;
        weights_current = weights_new;
        nll_previous    = nll_new;

        if outer >= 2 && max_dw <= weight_tol && abs(nll_change) <= nll_tol
            converged = true;
            break;
        end
        if no_improve_count >= no_improve_patience
            stopped_no_improve = true;
            break;
        end
    end

    if converged
        termination_reason = "converged_tolerances";
    elseif stopped_no_improve
        termination_reason = "stopped_no_likelihood_improvement";
    else
        termination_reason = "maximum_outer_iterations";
    end
end


function [Ab, taub, tsf_vector, tsf_matrix, fixed_nll, selected_nll, ...
    best_iteration, evaluated_iterations, termination_reason, error_message] = ...
    run_bootstrap_replicate(bootstrap_index, random_seed, chron_data, ...
        ok_idx, pz, Nz, pair_lo, pair_hi, pair_gap, w_order, delta_min_Ma, ...
        theta_fixed_hat, opts_boot, tsf_mode, group_indices, group_fixed, ...
        age_grid_step, hypsometry_pull, smoothness, update_fraction, ...
        smooth_span, max_outer, weight_tol, nll_tol, no_improve_patience)
% Run one independently seeded joint bootstrap resample.
    n_chron = numel(chron_data);
    Ab = nan(Nz, n_chron);
    taub = nan(n_chron, 1);
    tsf_vector = nan(Nz, 1);
    tsf_matrix = nan(Nz, n_chron);
    fixed_nll = nan;
    selected_nll = nan;
    best_iteration = nan;
    evaluated_iterations = nan;
    termination_reason = "";
    error_message = "";

    try
        % A dedicated substream makes replicate b identical in serial and
        % parallel execution, independent of worker scheduling.
        stream_b = RandStream('Threefry', 'Seed', random_seed);
        stream_b.Substream = bootstrap_index;

        cd_boot = chron_data;
        for c = ok_idx
            N_c = chron_data(c).N;
            idx_b = randi(stream_b, N_c, N_c, 1);
            cd_boot(c).age = chron_data(c).age(idx_b);
            cd_boot(c).sig = chron_data(c).sig(idx_b);
        end

        obj_fixed_b = @(th) nll_joint(th, cd_boot, ok_idx, pz, Nz, ...
            pair_lo, pair_hi, pair_gap, w_order, delta_min_Ma);
        [theta_b, fixed_nll] = fminunc( ...
            obj_fixed_b, theta_fixed_hat, opts_boot);
        selected_nll = fixed_nll;

        if tsf_mode == "iterative"
            [theta_b, tsf_matrix, selected_nll, best_iteration, ...
                evaluated_iterations, termination_reason] = ...
                fit_iterative_source_weights_bootstrap( ...
                    theta_b, fixed_nll, cd_boot, ok_idx, pz, Nz, ...
                    pair_lo, pair_hi, pair_gap, w_order, delta_min_Ma, ...
                    opts_boot, group_indices, group_fixed, age_grid_step, ...
                    hypsometry_pull, smoothness, update_fraction, ...
                    smooth_span, max_outer, weight_tol, nll_tol, ...
                    no_improve_patience);
        else
            tsf_vector = pz;
            best_iteration = 0;
            evaluated_iterations = 0;
            termination_reason = "fixed_source_weighting";
        end

        [Ab, taub] = unpack_joint(theta_b, cd_boot, ok_idx, Nz);
    catch ME
        error_message = string(ME.message);
    end
end


function report_bootstrap_progress(n_boot, reset_counter)
% Report ordered completion counts for serial or parallel bootstrap runs.
    persistent n_complete start_time
    if reset_counter || isempty(n_complete)
        n_complete = 0;
        start_time = tic;
        return;
    end

    n_complete = n_complete + 1;
    report_every = max(1, ceil(n_boot / 20));
    if n_boot <= 20 || mod(n_complete, report_every) == 0 || ...
            n_complete == n_boot
        elapsed_min = toc(start_time) / 60;
        remaining_min = elapsed_min / n_complete * (n_boot - n_complete);
        fprintf("  Bootstrap %d / %d  |  elapsed %.1f min  |  estimated remaining %.1f min\n", ...
            n_complete, n_boot, elapsed_min, remaining_min);
    end
end


function lookup = build_tc_lookup(TC)
% Build cell array of {name_variants, Tc_value} rows.
% Add more rows here to support additional chronometers.
    lookup = { ...
        {'ahe','aphe','ap-he','apatitehe','apatite-he','apatite_he', ...
         'apatite(u-th)he','apatite_u_th_he'},          TC.ahe; ...
        {'zhe','zirc-he','zirconhe','zircon-he','zircon_he', ...
         'zircon(u-th)he','zircon_u_th_he'},             TC.zhe; ...
        {'apb','appb','ap-pb','apatitepb','apatite-pb','apatite_pb', ...
         'apatite(u-th)pb','apatiteu-thpb','apupb'},     TC.apb; ...
        {'har','hbl','hornblende','hornblendear', ...
         'hornblende-ar','hornblende_ar','hblar', ...
         'hornblende40ar39ar'},                           TC.har  ...
    };
end
