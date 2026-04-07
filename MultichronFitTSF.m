%% MultichronFitTSF.m
%
% Joint multi-chronometer inversion of detrital age distributions using a
% hypsometry-weighted Gaussian mixture likelihood. All chronometers are
% optimized simultaneously with a closure-temperature-aware ordering penalty
% that enforces the physical age hierarchy (lower-Tc systems yield younger
% ages) without requiring explicit kinetic models.
%
% Supersedes fit_multichron_tsf.m (v4, independent solver).
% Same folder structure, same config format, extended outputs.
%
% KEY DIFFERENCE FROM v4
%   v4 fits each chronometer independently with fminsearch.
%   This script fits all chronometers jointly with fminunc (quasi-Newton),
%   adding a pairwise ordering penalty across all chronometer pairs scaled
%   by their normalized closure temperature separation. This prevents
%   physically implausible age inversions (e.g. ZHe older than APb) while
%   remaining agnostic about the specific thermal history.
%
% USAGE
%   1. Set catchment_name and base_dir below (the only two lines you edit).
%   2. Optionally override w_order / delta_min_Ma in the config CSV.
%   3. Run. All outputs are written into the catchment subfolder.
%
% FOLDER STRUCTURE  (unchanged from v4)
%   <base_dir>/
%   ├── MultichronFitTSF.m
%   ├── WP/
%   │   ├── WP_config.csv
%   │   ├── WP_Hypso.csv
%   │   ├── WP_ApHe.csv  (columns: Date_Ma, Error_Ma)
%   │   ├── WP_ZHe.csv
%   │   ├── WP_ApPb.csv
%   │   ├── WP_Hbl.csv
%   │   └── figures_svg/
%   └── ...
%
% CONFIG CSV FORMAT  (extended from v4)
%   Required columns (unchanged):
%     Chronometer, File, TauMin, AgeMargin, Lambda, AgeMinFilter, AgeMaxFilter
%   Optional global-override columns (add as single-row header + value, or
%   place in the Hypsometry row under these column names):
%     w_order     : ordering penalty weight  (default: see ORDERING PENALTY below)
%     delta_min   : min age separation scale (default: see ORDERING PENALTY below)
%
%   Example:
%     Chronometer,File,TauMin,AgeMargin,Lambda,AgeMinFilter,AgeMaxFilter
%     Hypsometry,WP_Hypso.csv,,,,,
%     ApHe,WP_ApHe.csv,0.5,5,0.0,0,Inf
%     ZHe,WP_ZHe.csv,0.5,5,0.0,0,Inf
%     ApPb,WP_ApPb.csv,0.1,15,0.5,0,Inf
%     Hbl,WP_Hbl.csv,0.1,5,1.0,83,Inf
%
% CLOSURE TEMPERATURES (Tc)
%   Hardcoded defaults from Hodges (2014) Table 2.
%   Edit the TC_DEFAULTS struct below to override for non-standard systems.
%   Recognised name variants are listed in build_tc_lookup().
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
%   predicted_bedrock_transect_<Chron>.csv          (same format as v4)
%   predicted_bedrock_transect_<Chron>_CI.csv       (bootstrap CI)
%   grain_posteriors_<Chron>.csv                    (optional)
%   grain_expected_source_<Chron>.csv
%   ordering_violations.csv                         (NEW: pre-fit violations)
%   ordering_penalty_contributions.csv              (NEW: per-pair penalty)
%   summary_fit_params.csv                          (extended from v4)
%   figures_svg/
%
% CHANGELOG
%   v5 : Joint fminunc optimizer with Tc-scaled pairwise ordering penalty.
%        Name-variant lookup for flexible chronometer labelling.
%        Ordering diagnostics (violation log, penalty contributions).
%        w_order / delta_min exposed as tunable parameters with config override.
%        Bootstrap calls joint solver per resample.
%   v4 : Config-file-driven workflow. Per-chronometer age filters.
%   v3 : Per-chronometer settings struct. Summary CSV.
%   v2 : Bug fixes (interp1, zu shadowing, cummax, tau init, quantile).

clear; close all; clc;

%% ============================================================
%% USER SETTINGS  <-- only edit these two lines between runs
% ---- Match these two lines to your MultichronFitTSF.m settings ----
catchment_name = "EX";                          % <-- change to your catchment name
base_dir       = "/path/to/your/project/folder"; % <-- change to your project folder
%% ============================================================


%% ---- CLOSURE TEMPERATURES  <-- edit here to override defaults ----
%
% Values from Hodges (2014) Table 2, bulk closure temperature (Tcb).
% Units: degrees Celsius.
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


%% ---- Global fitting options (rarely need changing) ----
age_grid_step          = 0.25;   % PDF plot resolution (Ma)
write_grain_posteriors = true;   % set false to skip large grain posterior CSVs
show_optimizer_iters   = false;  % show fminunc iterations in console
do_bootstrap           = true;
n_boot                 = 200;    % bootstrap resamples (500 for final runs)
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
                 'Lambda','AgeMinFilter','AgeMaxFilter'};
for rc = required_cols
    if ~any(strcmpi(cfg.Properties.VariableNames, rc{1}))
        error("Config file is missing required column: %s", rc{1});
    end
end

fprintf("Loaded config : %s\n",   config_file);
fprintf("Catchment     : %s\n",   catchment_name);
fprintf("Catchment dir : %s\n\n", catchment_dir);

% ---- Check for optional config overrides on Hypsometry row ----
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

hypsometry_file = fullfile(catchment_dir, cfg.File(hyps_mask));
chron_cfg       = cfg(chron_mask, :);
n_chron         = height(chron_cfg);

fprintf("Hypsometry file : %s\n", hypsometry_file);
fprintf("Chronometers    : %s\n", strjoin(chron_cfg.Chronometer, ', '));
fprintf("w_order         : %.2f\n", w_order);
fprintf("delta_min_Ma    : %.2f\n\n", delta_min_Ma);


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
for c = 1:n_chron
    cname = lower(strtrim(char(chron_cfg.Chronometer(c))));
    for row = 1:size(tc_lookup, 1)
        if any(strcmpi(tc_lookup{row,1}, cname))
            Tc_vec(c) = tc_lookup{row,2};
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
        fprintf("  %-12s  Tc = %d degC\n", chron_cfg.Chronometer(c), Tc_vec(c));
    else
        fprintf("  %-12s  Tc = unrecognised (no ordering constraint)\n", ...
            chron_cfg.Chronometer(c));
    end
end

% Build pairwise gap matrix: gap(i,j) = (Tc_j - Tc_i) / (Tc_max - Tc_min)
% Only populated where i < j and Tc_i < Tc_j (lower-Tc must be younger).
Tc_known  = Tc_vec(has_tc);
Tc_range  = max(Tc_known) - min(Tc_known);
if Tc_range == 0; Tc_range = 1; end  % guard single-chronometer edge case

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
        if Tc_vec(ci) < Tc_vec(cj)
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
    fprintf("  %s (Tc=%d) < %s (Tc=%d)  |  gap=%.3f  |  min_sep=%.3f Ma\n", ...
        chron_cfg.Chronometer(pair_lo(p)), Tc_vec(pair_lo(p)), ...
        chron_cfg.Chronometer(pair_hi(p)), Tc_vec(pair_hi(p)), ...
        pair_gap(p), delta_min_Ma * pair_gap(p));
end
fprintf("\n");


%% ---------------- READ ALL DETRITAL DATA ----------------
% Load all chronometer data first so the joint objective can access them.
fprintf("Loading detrital data...\n");
chron_data = struct();   % chron_data(c).age, .sig, .N, .AminB, .AmaxB, etc.

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

fprintf("Joint fit complete.  Total NLL = %.4f\n\n", nll_total);

% Unpack joint solution
[Ahat_all, tauhat_all] = unpack_joint(theta_joint_hat, chron_data, ok_idx, Nz);


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

Ab_all   = nan(Nz, n_chron, n_boot);
taub_all = nan(n_chron, n_boot);

if do_bootstrap
    fprintf("Phase 3: joint bootstrap (%d resamples)...\n", n_boot);
    fprintf("  (Each resample is a full joint fminunc solve -- may take a few minutes)\n");

    opts_boot = optimoptions('fminunc', ...
        'Algorithm',              'quasi-newton', ...
        'MaxIterations',          3000, ...
        'MaxFunctionEvaluations', 5e5, ...
        'FiniteDifferenceType',   'central', ...
        'Display',                'off', ...
        'OptimalityTolerance',    1e-7, ...
        'StepTolerance',          1e-9);

    % Warm-start from joint solution
    theta_warm = theta_joint_hat;

    for b = 1:n_boot
        % Resample each chronometer independently, build resampled data struct
        cd_boot = chron_data;
        for c = ok_idx
            N_c       = chron_data(c).N;
            idx_b     = randi(N_c, N_c, 1);
            cd_boot(c).age = chron_data(c).age(idx_b);
            cd_boot(c).sig = chron_data(c).sig(idx_b);
        end

        obj_b = @(th) nll_joint(th, cd_boot, ok_idx, pz, Nz, ...
                                 pair_lo, pair_hi, pair_gap, ...
                                 w_order, delta_min_Ma);

        try
            th_b = fminunc(obj_b, theta_warm, opts_boot);
            [Ab, taub] = unpack_joint(th_b, chron_data, ok_idx, Nz);
            Ab_all(:,:,b)  = Ab;
            taub_all(:,b)  = taub;
        catch
            % Resample failed; leave as NaN (excluded from quantiles)
        end

        if mod(b, 50) == 0
            fprintf("  Bootstrap %d / %d\n", b, n_boot);
        end
    end
    fprintf("  Bootstrap complete.\n\n");
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

    fprintf("===== %s  (joint fit) =====\n", chron);
    fprintf("  tau = %.3f Ma  |  A = [%.2f, %.2f] Ma\n", ...
        tauhat, min(Ahat), max(Ahat));

    % NLL for this chronometer alone (for summary reporting)
    nll_c = nll_single(theta_joint_hat( block_slice(c, ok_idx, Nz) ), ...
                       age, sig, pz, cd.AminB, cd.AmaxB, cd.tau_min, cd.lambda);
    fprintf("  Single-chron NLL contribution = %.4f\n\n", nll_c);

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

        out_ci = table(z_centers, pz, Ahat, A_med, A_lo, A_hi, ...
            A_sigma_ish, A_sigma_minus, A_sigma_plus, ...
            'VariableNames', {'Elevation_m','HypsometryWeight','Ahat_Ma', ...
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
        y_pred = y_pred + pz(k) * normpdf(age_grid, Ahat(k), sig_pred);
    end
    y_pred = y_pred / trapz(age_grid, y_pred);

    figure('Name', sprintf('%s_PDF_Obs_vs_Pred', chron)); hold on;
    plot(age_grid, y_obs, 'LineWidth', 1.5);
    plot(age_grid, y_pred, 'LineWidth', 1.5);
    xlabel('Age (Ma)'); ylabel('Probability density');
    legend('Observed (kernel)', 'Predicted (hypsometry mixture)', 'Location', 'best');
    title(sprintf('%s: observed vs predicted detrital PDF', chron));

    %% ---- Per-grain posterior p(z_k | age_i) ----
    Pzik = zeros(Nz, N);
    for i = 1:N
        s      = sqrt(sig(i)^2 + tauhat^2);
        like_k = normpdf(age(i), Ahat, s);
        post   = like_k(:) .* pz(:);
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
    title(sprintf('%s: joint-fit age-elevation curve', chron));

    figure('Name', sprintf('%s_CDF_Hyps_vs_Implied', chron)); hold on;
    stairs(z_edges(2:end), cumsum(pz),         'LineWidth', 1.5);
    stairs(z_edges(2:end), cumsum(pz_implied), 'LineWidth', 1.5);
    xlabel('Elevation (m)'); ylabel('Cumulative probability');
    legend('Hypsometry CDF (TSF)', 'Implied source CDF', 'Location', 'best');
    title(sprintf('%s: hypsometry vs implied source distribution', chron));

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
    out_transect = table(z_centers, pz, Ahat, ...
        'VariableNames', {'Elevation_m','HypsometryWeight','PredictedBedrockAge_Ma'});
    writetable(out_transect, fullfile(catchment_dir, ...
        "predicted_bedrock_transect_" + string(chron) + ".csv"));
    fprintf("  Wrote transect CSV\n");

    if write_grain_posteriors
        Pfile  = fullfile(catchment_dir, "grain_posteriors_" + string(chron) + ".csv");
        Ptable = array2table(Pzik);
        Ptable = addvars(Ptable, z_centers, 'Before', 1, 'NewVariableNames', "Elevation_m");
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
    Summary = [Summary; table(chron_str, N, cd.n_excluded, nll_c, tauhat, ...
        min(Ahat), max(Ahat), cd.tau_min, cd.age_margin, cd.lambda, ...
        cd.age_min_filt, cd.age_max_filt, w_order, delta_min_Ma, Tc_vec(c), ...
        'VariableNames', {'Chronometer','Ngrains','Nexcluded','NLL_single', ...
                          'Tau_Ma','Amin_Ma','Amax_Ma','Setting_TauMin', ...
                          'Setting_AgeMargin','Setting_Lambda', ...
                          'Setting_AgeMinFilter','Setting_AgeMaxFilter', ...
                          'Setting_w_order','Setting_delta_min_Ma','Tc_degC'})]; %#ok<AGROW>

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
    title(sprintf('%s: all chronometers — joint fit', catchment_name));
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
% NLL for a single chronometer (unchanged from v4 except name).
    Nz = numel(pz);
    [A, tau] = unpack_single(theta, Nz, AminB, AmaxB, tau_min);
    nll = 0;
    for i = 1:numel(age)
        s    = sqrt(sig(i)^2 + tau^2);
        comp = pz(:) .* normpdf(age(i), A(:), s);
        p    = max(sum(comp), realmin);
        nll  = nll - log(p);
    end
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


function nll = nll_joint(theta, chron_data, ok_idx, pz, Nz, ...
                          pair_lo, pair_hi, pair_gap, w_order, delta_min_Ma)
% Joint NLL: sum of per-chronometer NLLs + pairwise ordering penalty.

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

        nll = nll + nll_single(th_c, cd.age, cd.sig, pz, ...
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
