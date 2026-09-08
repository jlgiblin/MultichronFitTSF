%% MultichronFitTSF_Georef.m
%
% DEM-BASED COORDINATE ASSIGNMENT MODULE
% =========================================
% Reads a GeoTIFF DEM and a flow accumulation GeoTIFF, identifies stream/channel
% pixels within each hypsometric elevation bin, and assigns representative
% lat/lon coordinates (projected to the channel pixel closest to the median
% elevation of each bin) to your existing predicted bedrock transect outputs.
%
% Designed to slot directly into the MultichronFitTSF workflow.
% Run AFTER MultichronFitTSF.m has produced its transect CSV outputs.
%
% ---- REQUIRED INPUTS ----
%   1. DEM GeoTIFF          : Elevation raster clipped to catchment (.tif)
%   2. Flow accumulation    : GeoTIFF of upstream contributing area (pixels),
%                             derived from DEM in ArcGIS/QGIS/MATLAB.
%                             Pixels with value > flow_acc_threshold are "channels".
%   3. Hypsometry CSV       : Same file used in the fitting script
%                             (columns: Elevation_m or Elevation, RelArea or Area)
%   4. Transect CSVs        : Output CSVs from MultichronFitTSF.m
%                             e.g. predicted_bedrock_transect_ApHe_CI.csv
%                             Must have column: Elevation_m
%   5. Grain source CSVs    : Output from MultichronFitTSF.m
%                             e.g. grain_expected_source_ApHe.csv
%                             Must have columns: Date_Ma, Median_SourceElev_m, etc.
%
% ---- OUTPUTS ----
%   - predicted_bedrock_transect_<Chron>_georef.csv
%       Same as input transect CSV plus the requested geographic and/or
%       projected coordinates, Channel_Elev_m, and N_channel_pixels_in_bin
%   - grain_expected_source_<Chron>_georef.csv
%       Same as grain source CSV plus median and P16/P84 source coordinates
%   - elevation_coordinate_assignment.csv
%       One row per elevation bin with coordinate/elevation mismatch QA fields
%   - <catchment>_dem_coord_assignment.svg/pdf/png
%       Map showing DEM hillshade, channel network, and numbered representatives
%
% ---- TOOLBOX REQUIREMENTS ----
%   MATLAB Mapping Toolbox (for readgeoraster, geotiffinfo, projfwd/projinv)
%   OR Image Processing Toolbox as fallback (geotiffread - older API)
%   No additional toolboxes required for core computation.
%
% ---- HOW TO DERIVE FLOW ACCUMULATION ----
%   ArcGIS:  Spatial Analyst > Hydrology > Flow Direction, then Flow Accumulation
%   QGIS:    SAGA > Terrain Analysis > Flow Accumulation (top-down)
%   MATLAB:  See note at bottom of this file using TopoToolbox (optional)
%
% 

clear; close all; clc;

%% ============================================================
%% USER SETTINGS  <-- edit this block only
%% ============================================================

% ---- Match these two lines to your MultichronFitTSF.m settings ----
catchment_name = "EX";                          % <-- change to your catchment name
base_dir       = "/path/to/your/project/folder"; % <-- change to your project folder

% ---- DEM and flow accumulation GeoTIFFs ----
% By default, assumed to live in the same catchment subfolder as your data.
% Override with full paths if stored elsewhere.
catchment_dir = fullfile(base_dir, catchment_name);
config_file   = fullfile(catchment_dir, catchment_name + "_config.csv");
dem_file      = fullfile(catchment_dir, catchment_name + "_DEM.tif");      % clipped DEM GeoTIFF
flowacc_file  = fullfile(catchment_dir, catchment_name + "_flowacc.tif");  % flow accumulation GeoTIFF

% ---- Channel definition ----
% Pixels with flow accumulation > this value are treated as channels.
% There is no universal value: increase it to retain fewer, larger channels;
% decrease it to include more headwater pixels. Inspect the diagnostic map
% and elevation_coordinate_assignment.csv after every change.
flow_acc_threshold = 100;

% Warn when the assigned channel pixel differs from its modeled elevation by
% more than this amount. This is a QA limit, not an alternate search rule.
max_elevation_mismatch_m = 50;

% ---- Coordinate output format ----
% "geographic"  => output Lat/Lon (WGS84 decimal degrees)
% "projected"   => output Easting/Northing in the DEM's native CRS
% "both"        => output all four columns
coord_output_mode = "both";

% ---- Elevation bin matching tolerance ----
% When assigning a channel pixel to an elevation bin, allow pixels within
% +/- this many meters of the modeled bin center to count. The highest-flow
% pixel in that narrow elevation window is selected. If none exists, the
% nearest channel pixel is used and explicitly flagged for review.
elev_tolerance_m = 10;

% ---- Output directory ----
% Outputs are written to a georef_outputs/ subfolder inside the catchment folder.
% Override with a full path if you prefer a different location.
output_dir = fullfile(catchment_dir, "georef_outputs");

%% ============================================================
%% END USER SETTINGS
%% ============================================================

% Read the fitting config so chronometer labels and the hypsometry filename
% always match the MultichronFitTSF run. Users do not maintain a second list.
if ~isfile(config_file)
    error("Config file not found: %s", config_file);
end
cfg = readtable(config_file, 'TextType', 'string');
if ~all(ismember({'Chronometer','File'}, cfg.Properties.VariableNames))
    error("Config file must contain Chronometer and File columns: %s", config_file);
end
hyps_mask = strcmpi(strtrim(cfg.Chronometer), 'Hypsometry');
if nnz(hyps_mask) ~= 1
    error("Config file must contain exactly one Hypsometry row: %s", config_file);
end
hypsometry_csv = fullfile(catchment_dir, cfg.File(hyps_mask));

chron_names = strtrim(cfg.Chronometer(~hyps_mask));
chron_names = chron_names(~ismissing(chron_names) & strlength(chron_names) > 0);
if isempty(chron_names)
    error("No chronometer rows were found in: %s", config_file);
end

% Prefer bootstrap-CI transects when present, but also support preliminary
% non-bootstrap runs. Grain source files use the same config-derived labels.
transect_files = strings(size(chron_names));
for c = 1:numel(chron_names)
    ci_candidate = fullfile(catchment_dir, ...
        "predicted_bedrock_transect_" + chron_names(c) + "_CI.csv");
    fit_candidate = fullfile(catchment_dir, ...
        "predicted_bedrock_transect_" + chron_names(c) + ".csv");
    if isfile(ci_candidate)
        transect_files(c) = ci_candidate;
    else
        transect_files(c) = fit_candidate;
    end
end
grain_files = fullfile(catchment_dir, ...
    "grain_expected_source_" + chron_names + ".csv");

if ~exist(output_dir, "dir"), mkdir(output_dir); end

fprintf("=== DEM Coordinate Assignment Module ===\n");
fprintf("Catchment: %s\n\n", catchment_name);

%% ---- STEP 1: Read DEM and flow accumulation rasters ----
fprintf("Reading DEM: %s\n", dem_file);
[dem_Z, dem_R] = readDEM_safe(dem_file);

fprintf("Reading flow accumulation: %s\n", flowacc_file);
[acc_Z, acc_R] = readDEM_safe(flowacc_file);

% Sanity check: rasters should be same size
if ~isequal(size(dem_Z), size(acc_Z))
    error("DEM and flow accumulation rasters have different sizes: " + ...
          "DEM is %dx%d, FlowAcc is %dx%d. " + ...
          "Please ensure both are clipped to the same extent and resolution.", ...
           size(dem_Z,1), size(dem_Z,2), size(acc_Z,1), size(acc_Z,2));
end
if ~rasterReferencesAlign(dem_R, acc_R)
    error("DEM and flow accumulation rasters do not share the same spatial " + ...
          "reference. Reproject, resample, and clip them to an identical grid.");
end

fprintf("DEM size: %d rows x %d cols\n", size(dem_Z,1), size(dem_Z,2));
fprintf("DEM elevation range: %.1f – %.1f m\n", ...
    min(dem_Z(:), [],'omitnan'), max(dem_Z(:), [],'omitnan'));

%% ---- STEP 2: Build channel mask ----
channel_mask = (acc_Z > flow_acc_threshold) & ~isnan(dem_Z);
n_channel_px = sum(channel_mask(:));
fprintf("Channel pixels (flow acc > %d): %d of %d total\n", ...
    flow_acc_threshold, n_channel_px, sum(~isnan(dem_Z(:))));

if n_channel_px == 0
    error("No channel pixels found with flow_acc_threshold = %d. " + ...
          "Lower the threshold or check that your flow accumulation " + ...
          "raster covers the same area as the DEM.", flow_acc_threshold);
end

%% ---- STEP 3: Get geographic coordinates of all DEM pixels ----
fprintf("Computing pixel coordinates...\n");
[nrows, ncols] = size(dem_Z);

% Build row/col grids
[col_grid, row_grid] = meshgrid(1:ncols, 1:nrows);

% Convert pixel row/col -> geographic coordinates using the raster reference
[lat_grid, lon_grid] = pixel2latlon(dem_R, row_grid, col_grid);

% Also get projected coordinates when the DEM uses a projected reference.
is_projected_dem = isa(dem_R, 'map.rasterref.MapCellsReference') || ...
                   isa(dem_R, 'map.rasterref.MapPostingsReference');
if is_projected_dem
    [x_grid, y_grid] = pixel2xy(dem_R, row_grid, col_grid);
else
    x_grid = nan(size(lat_grid));
    y_grid = nan(size(lat_grid));
end

coord_output_mode = lower(strtrim(coord_output_mode));
if ~any(coord_output_mode == ["geographic", "projected", "both"])
    error("coord_output_mode must be geographic, projected, or both.");
end
if coord_output_mode == "projected" && ~is_projected_dem
    error("coord_output_mode is projected, but the DEM has a geographic " + ...
          "reference. Use a projected DEM or select geographic output.");
end
include_geographic = coord_output_mode == "geographic" || coord_output_mode == "both";
include_projected  = (coord_output_mode == "projected" || ...
                      coord_output_mode == "both") && is_projected_dem;
if coord_output_mode == "both" && ~is_projected_dem
    warning("The DEM is geographic, so projected Easting/Northing columns " + ...
            "cannot be calculated. Writing Lat/Lon only.");
end

%% ---- STEP 4: Read hypsometry and reconstruct elevation bins ----
fprintf("Reading hypsometry: %s\n", hypsometry_csv);
H = readtable(hypsometry_csv);

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
    A_raw(A_raw < 0) = 0;
    F_raw = A_raw / A_raw(end);
else
    error("Hypsometry file must contain 'RelArea' or 'Area'.");
end

[z_raw, idx] = sort(z_raw);
F_raw = F_raw(idx);
[zu, ia] = unique(z_raw, 'stable');
Fu = F_raw(ia);
Fu = (Fu - Fu(1)) / (Fu(end) - Fu(1));
Fu(1) = 0; Fu(end) = 1;

% Infer the number of equal-area bins from an actual transect output. This
% remains correct if target_bins was changed in MultichronFitTSF.m.
reference_transect_idx = find(isfile(transect_files), 1);
if isempty(reference_transect_idx)
    error("No predicted_bedrock_transect_<Chron>[_CI].csv files were found in %s", ...
        catchment_dir);
end
reference_transect = readtable(transect_files(reference_transect_idx));
if ~any(strcmpi(reference_transect.Properties.VariableNames, 'Elevation_m'))
    error("Reference transect is missing Elevation_m: %s", ...
        transect_files(reference_transect_idx));
end
target_bins = height(reference_transect);
if target_bins < 2
    error("Reference transect must contain at least two elevation bins: %s", ...
        transect_files(reference_transect_idx));
end

Fq = linspace(0, 1, target_bins+1)';
[Fu_u, iu] = unique(Fu, 'stable');
zu_u = zu(iu);
z_edges = interp1(Fu_u, zu_u, Fq, 'linear', 'extrap');
z_edges = max(z_edges, cummax(z_edges));
z_edges(1) = zu_u(1); z_edges(end) = zu_u(end);

z_centers = 0.5 * (z_edges(1:end-1) + z_edges(2:end));
Nz = numel(z_centers);

reference_centers = reference_transect.Elevation_m(:);
center_tolerance = max(1e-6, 1e-8 * max(abs(reference_centers)));
if numel(reference_centers) ~= Nz || ...
        max(abs(reference_centers - z_centers)) > center_tolerance
    error("Hypsometry and transect elevation bins do not match. Confirm that " + ...
          "the hypsometry file is the one used to create the transect outputs.");
end
z_centers = reference_centers;

fprintf("Reconstructed %d elevation bins: %.0f – %.0f m\n", ...
    Nz, z_edges(1), z_edges(end));

%% ---- STEP 5: Assign channel pixel representative to each elevation bin ----
fprintf("Assigning channel pixels to elevation bins...\n");

bin_lat      = nan(Nz, 1);
bin_lon      = nan(Nz, 1);
bin_easting  = nan(Nz, 1);
bin_northing = nan(Nz, 1);
bin_channel_elev = nan(Nz, 1);
bin_n_channel_px = zeros(Nz, 1);
bin_used_fallback = false(Nz, 1);
bin_flow_acc = nan(Nz, 1);

% Flatten grids for channel pixels
chan_idx  = find(channel_mask);
chan_elev = dem_Z(chan_idx);
chan_lat  = lat_grid(chan_idx);
chan_lon  = lon_grid(chan_idx);
chan_x    = x_grid(chan_idx);
chan_y    = y_grid(chan_idx);
chan_acc  = acc_Z(chan_idx);   % flow accumulation value

for k = 1:Nz
    in_bin = abs(chan_elev - z_centers(k)) <= elev_tolerance_m;

    if ~any(in_bin)
        % Widen search: find closest channel pixel by elevation
        [~, closest] = min(abs(chan_elev - z_centers(k)));
        in_bin(closest) = true;
        bin_used_fallback(k) = true;
        fprintf("  Bin %d (%.0f m): no channel pixel in range — using nearest (%.0f m away)\n", ...
            k, z_centers(k), abs(chan_elev(closest) - z_centers(k)));
    end

    bin_n_channel_px(k) = sum(in_bin);

    % Within the narrow elevation window, prefer the pixel with the highest
    % flow accumulation to favor the primary sediment-routing path.
    in_bin_idx = find(in_bin);
    [~, best] = max(chan_acc(in_bin_idx));
    rep_idx = in_bin_idx(best);

    bin_lat(k)          = chan_lat(rep_idx);
    bin_lon(k)          = chan_lon(rep_idx);
    bin_easting(k)      = chan_x(rep_idx);
    bin_northing(k)     = chan_y(rep_idx);
    bin_channel_elev(k) = chan_elev(rep_idx);
    bin_flow_acc(k)     = chan_acc(rep_idx);
end

fprintf("Coordinate assignment complete.\n\n");

bin_elev_difference = abs(bin_channel_elev - z_centers);
bin_exceeds_mismatch = bin_elev_difference > max_elevation_mismatch_m;

% Write one compact QA table before repeating these coordinates in each
% chronometer-specific transect.
Assignment = table((1:Nz)', z_centers, bin_channel_elev, ...
    bin_elev_difference, bin_flow_acc, bin_n_channel_px, ...
    bin_used_fallback, bin_exceeds_mismatch, ...
    'VariableNames', {'Bin','Elevation_m','Channel_Elev_m', ...
                      'Elevation_Difference_m','FlowAccumulation', ...
                      'N_channel_px','UsedFallback','ExceedsMismatchLimit'});
if include_geographic
    Assignment.Lat = bin_lat;
    Assignment.Lon = bin_lon;
end
if include_projected
    Assignment.Easting_m = bin_easting;
    Assignment.Northing_m = bin_northing;
end
assignment_file = fullfile(output_dir, "elevation_coordinate_assignment.csv");
writetable(Assignment, assignment_file);
fprintf("Wrote coordinate-assignment QA table: %s\n", assignment_file);

if any(bin_exceeds_mismatch)
    warning("%d elevation bins exceed the %.1f m coordinate-assignment " + ...
            "mismatch limit. Lower flow_acc_threshold, inspect the " + ...
            "diagnostic map, or verify the raster coverage.", ...
        nnz(bin_exceeds_mismatch), max_elevation_mismatch_m);
end

% Report any bins with very few channel pixels (may indicate DEM/bin mismatch)
sparse_bins = find(bin_n_channel_px < 3);
if ~isempty(sparse_bins)
    fprintf("WARNING: %d bins have fewer than 3 channel pixels:\n", numel(sparse_bins));
    for kb = sparse_bins(:)'
        fprintf("  Bin %d: center=%.0f m, n_channel=%d, assigned lat=%.5f lon=%.5f\n", ...
            kb, z_centers(kb), bin_n_channel_px(kb), bin_lat(kb), bin_lon(kb));
    end
    fprintf("Consider lowering flow_acc_threshold or checking DEM coverage.\n\n");
end

%% ---- STEP 6: Save diagnostic figure ----
fprintf("Generating diagnostic map...\n");

% --- Hillshade from DEM gradient ---
[dzdx, dzdy] = gradient(dem_Z);
hillshade = -dzdx * cosd(315)*cosd(45) - dzdy * sind(315)*cosd(45) + sind(45);
hillshade(isnan(dem_Z)) = NaN;

% --- Normalize hillshade and elevation to [0,1] ---
hs_min = min(hillshade(:), [], 'omitnan'); hs_max = max(hillshade(:), [], 'omitnan');
el_min = min(dem_Z(:),    [], 'omitnan'); el_max = max(dem_Z(:),    [], 'omitnan');
hs_norm = (hillshade - hs_min) / (hs_max - hs_min);
el_norm = (dem_Z    - el_min)  / (el_max - el_min);

% --- Weighted blend ---
hs_weight = 0.45; el_weight = 0.55;
blend = hs_weight * hs_norm + el_weight * el_norm;

% --- Fill NaN (outside catchment) with a flat value so no transparency is needed ---
% 0 will map to the first colormap entry, which we set to neutral gray below
blend_filled = blend;
blend_filled(isnan(dem_Z)) = 0;

% --- Shift catchment values to [0.05, 1] so gray (0) is visually distinct ---
blend_filled(~isnan(dem_Z)) = 0.05 + 0.95 * blend(~isnan(dem_Z));

% --- Custom colormap: index 1 = gray (outside catchment), rest = parula ---
n_colors = 256;
gray_entry = [0.75 0.75 0.75];
cmap = [gray_entry; parula(n_colors - 1)];

% --- Figure setup ---
fig_diag = figure('Name', 'DEM_coord_assignment', ...
                  'Units', 'inches', 'Position', [1 1 8 6.5]);

ax = axes('Units', 'normalized', 'Position', [0.08 0.08 0.78 0.86]);

% Single imagesc, no AlphaData — flat raster, no transparency
imagesc(ax, blend_filled, [0 1]);
colormap(ax, cmap);

% --- Colorbar with ticks remapped to real elevation values ---
cb = colorbar(ax, 'eastoutside');
cb.Label.String = 'Elevation (m)';
cb.Label.FontSize = 9;
% Ticks only over the catchment range [0.05, 1], not the gray zone
cb.Limits = [0.05 1];
cb.Ticks = linspace(0.05, 1, 6);
cb.TickLabels = arrayfun(@(v) sprintf('%.0f', ...
    el_min + ((v - 0.05) / 0.95) * (el_max - el_min)), ...
    cb.Ticks, 'UniformOutput', false);

set(ax, 'YDir', 'normal', 'FontSize', 8, 'TickDir', 'out', 'Box', 'off');
axis(ax, 'equal', 'tight');
hold(ax, 'on');

% --- Channel network ---
[cr, cc] = find(channel_mask);
plot(ax, cc, cr, '.', 'Color', [0.15 0.45 0.85], ...
     'MarkerSize', 1.5, 'DisplayName', 'Channel network');

% --- Bin representative points and labels ---
for k = 1:Nz
    if ~isnan(bin_lat(k))
        [pr, pc] = latlon2pixel(dem_R, bin_lat(k), bin_lon(k));
        scatter(ax, pc, pr, 55, 'filled', ...
                'MarkerFaceColor', 'w', ...
                'MarkerEdgeColor', [0.15 0.15 0.15], ...
                'LineWidth', 0.8, ...
                'HandleVisibility', 'off');
        text(ax, pc, pr, sprintf('%d', k), ...
             'FontSize', 6, 'FontWeight', 'bold', ...
             'Color', [0.05 0.05 0.05], ...
             'HorizontalAlignment', 'center', ...
             'VerticalAlignment', 'middle', ...
             'Interpreter', 'none');
    end
end

% --- Labels and legend ---
title(ax, sprintf('%s — elevation bin channel representatives', catchment_name), ...
      'FontSize', 10, 'FontWeight', 'normal', 'Interpreter', 'none');
xlabel(ax, 'Column (pixel)', 'FontSize', 8);
ylabel(ax, 'Row (pixel)',    'FontSize', 8);
legend(ax, 'Channel network', ...
       'Location', 'northeast', 'FontSize', 7, 'Box', 'off');

% --- Export ---
diag_fig_file_svg = fullfile(output_dir, catchment_name + "_dem_coord_assignment.svg");
print(fig_diag, diag_fig_file_svg, '-dsvg', '-vector');
fprintf("Saved vector diagnostic map: %s\n", diag_fig_file_svg);

diag_fig_file_pdf = fullfile(output_dir, catchment_name + "_dem_coord_assignment.pdf");
exportgraphics(fig_diag, diag_fig_file_pdf, 'ContentType', 'vector', 'Resolution', 300);
fprintf("Saved vector PDF diagnostic map: %s\n", diag_fig_file_pdf);

diag_fig_file_png = fullfile(output_dir, catchment_name + "_dem_coord_assignment.png");
exportgraphics(fig_diag, diag_fig_file_png, 'Resolution', 250);
fprintf("Saved raster diagnostic map: %s\n\n", diag_fig_file_png);
%% ---- STEP 7: Loop over chronometers — add coords to transect and grain files ----
for c = 1:numel(chron_names)
    chron = chron_names(c);
    fprintf("===== %s =====\n", chron);

    %% -- 7a: Update transect CSV --
    tfile = transect_files(c);
    if ~isfile(tfile)
        warning("Transect file not found: %s. Skipping.", tfile);
    else
        T = readtable(tfile);

        % Match each transect row to the nearest elevation bin
        if ~any(strcmpi(T.Properties.VariableNames, 'Elevation_m'))
            warning("Transect file %s missing 'Elevation_m' column. Skipping.", tfile);
        else
            T_elev = T.Elevation_m;
            [~, bin_idx] = min(abs(T_elev(:) - z_centers(:)'), [], 2);

            if include_geographic
                T.Lat = bin_lat(bin_idx);
                T.Lon = bin_lon(bin_idx);
            end
            if include_projected
                T.Easting_m  = bin_easting(bin_idx);
                T.Northing_m = bin_northing(bin_idx);
            end
            T.Channel_Elev_m = bin_channel_elev(bin_idx);
            T.Elevation_Difference_m = abs(T_elev(:) - T.Channel_Elev_m);
            T.FlowAccumulation = bin_flow_acc(bin_idx);
            T.N_channel_px  = bin_n_channel_px(bin_idx);
            T.UsedFallback  = bin_used_fallback(bin_idx);
            T.ExceedsMismatchLimit = bin_exceeds_mismatch(bin_idx);

            out_tfile = fullfile(output_dir, ...
                "predicted_bedrock_transect_" + chron + "_georef.csv");
            writetable(T, out_tfile);
            fprintf("  Wrote georeferenced transect: %s\n", out_tfile);
        end
    end

    %% -- 7b: Update grain source CSV --
    gfile = grain_files(c);
    if ~isfile(gfile)
        warning("Grain source file not found: %s. Skipping.", gfile);
    else
        G = readtable(gfile);

        % Assign coordinates based on the MEDIAN source elevation per grain
        % (using the Median_SourceElev_m column from your existing output)
        if ~any(strcmpi(G.Properties.VariableNames, 'Median_SourceElev_m'))
            warning("Grain file %s missing 'Median_SourceElev_m'. " + ...
                    "Trying 'E_SourceElev_m'.", gfile);
            if any(strcmpi(G.Properties.VariableNames, 'E_SourceElev_m'))
                grain_elev_col = G.E_SourceElev_m;
            else
                warning("Cannot find source elevation column in %s. Skipping.", gfile);
                continue;
            end
        else
            grain_elev_col = G.Median_SourceElev_m;
        end

        % Match grain elevation to nearest bin
        [~, grain_bin_idx] = min(abs(grain_elev_col(:) - z_centers(:)'), [], 2);

        if include_geographic
            G.Lat_median = bin_lat(grain_bin_idx);
            G.Lon_median = bin_lon(grain_bin_idx);
        end
        if include_projected
            G.Easting_median  = bin_easting(grain_bin_idx);
            G.Northing_median = bin_northing(grain_bin_idx);
        end

        % Also assign based on P16 and P84 source elevation for uncertainty range
        if all(ismember({'P16_SourceElev_m','P84_SourceElev_m'}, ...
                G.Properties.VariableNames))
            [~, g_p16] = min(abs(G.P16_SourceElev_m(:) - z_centers(:)'), [], 2);
            [~, g_p84] = min(abs(G.P84_SourceElev_m(:) - z_centers(:)'), [], 2);
            if include_geographic
                G.Lat_p16 = bin_lat(g_p16);
                G.Lon_p16 = bin_lon(g_p16);
                G.Lat_p84 = bin_lat(g_p84);
                G.Lon_p84 = bin_lon(g_p84);
            end
            if include_projected
                G.Easting_p16  = bin_easting(g_p16);
                G.Northing_p16 = bin_northing(g_p16);
                G.Easting_p84  = bin_easting(g_p84);
                G.Northing_p84 = bin_northing(g_p84);
            end
        end

        out_gfile = fullfile(output_dir, ...
            "grain_expected_source_" + chron + "_georef.csv");
        writetable(G, out_gfile);
        fprintf("  Wrote georeferenced grain source file: %s\n", out_gfile);
    end

    fprintf("\n");
end

fprintf("=== Done. All outputs written to: %s ===\n", output_dir);

%% ============================================================
%% LOCAL HELPER FUNCTIONS
%% ============================================================

function [Z, R] = readDEM_safe(filepath)
    % Reads a GeoTIFF using readgeoraster (Mapping Toolbox, R2020b+).
    % Falls back to geotiffread (older API) if needed.
    if ~isfile(filepath)
        error("File not found: %s", filepath);
    end
    try
        [Z, R] = readgeoraster(filepath);
        Z = double(Z);
        Z(Z < -1e10) = NaN;   % replace common NoData values
        Z(Z == 0 & all(Z(:)==0)) = NaN; % all-zero edge case
    catch ME1
        try
            warning("readgeoraster failed (%s). Trying geotiffread...", ME1.message);
            [Z, R] = geotiffread(filepath);
            Z = double(Z);
            Z(Z < -1e10) = NaN;
        catch ME2
            error("Cannot read GeoTIFF %s.\nreadgeoraster: %s\ngeotiffread: %s", ...
                filepath, ME1.message, ME2.message);
        end
    end
end

function tf = rasterReferencesAlign(R1, R2)
    % Confirm that two rasters occupy the same grid, not merely the same size.
    if ~strcmp(class(R1), class(R2)) || ~isequal(R1.RasterSize, R2.RasterSize)
        tf = false;
        return;
    end
    if isprop(R1, 'XWorldLimits') && isprop(R1, 'YWorldLimits')
        v1 = [R1.XWorldLimits, R1.YWorldLimits];
        v2 = [R2.XWorldLimits, R2.YWorldLimits];
    elseif isprop(R1, 'LatitudeLimits') && isprop(R1, 'LongitudeLimits')
        v1 = [R1.LatitudeLimits, R1.LongitudeLimits];
        v2 = [R2.LatitudeLimits, R2.LongitudeLimits];
    else
        tf = isequal(R1, R2);
        return;
    end
    scale = max(1, max(abs([v1, v2])));
    tf = max(abs(v1 - v2)) <= 1e-10 * scale;
end

function [lat, lon] = pixel2latlon(R, row, col)
    % Convert pixel row/col to geographic (lat/lon) using raster reference.
    % Handles both MapCellsReference (projected) and GeographicCellsReference.
    if isa(R, 'map.rasterref.MapCellsReference') || ...
       isa(R, 'map.rasterref.MapPostingsReference')
        % Projected CRS — get x/y then convert
        [x, y] = R.intrinsicToWorld(col, row);
        try
            proj = R.ProjectedCRS;
            [lat, lon] = projinv(proj, x, y);
        catch
            error("Projected DEM has no usable embedded CRS. Re-export the " + ...
                  "GeoTIFF with its CRS before requesting Lat/Lon output.");
        end
    elseif isa(R, 'map.rasterref.GeographicCellsReference') || ...
           isa(R, 'map.rasterref.GeographicPostingsReference')
        % Geographic CRS — direct conversion
        [lat, lon] = R.intrinsicToGeographic(col, row);
    else
        % Fallback: try intrinsicToWorld
        try
            [x, y] = R.intrinsicToWorld(col, row);
            lat = y; lon = x;
        catch
            error("Unrecognized raster reference type: %s", class(R));
        end
    end
end

function [x, y] = pixel2xy(R, row, col)
    % Get projected x/y (or geographic lon/lat) for pixel positions.
    try
        [x, y] = R.intrinsicToWorld(col, row);
    catch
        [y, x] = pixel2latlon(R, row, col);
    end
end

function [row, col] = latlon2pixel(R, lat, lon)
    % Convert lat/lon back to pixel row/col (for plotting only).
    if isa(R, 'map.rasterref.MapCellsReference') || ...
       isa(R, 'map.rasterref.MapPostingsReference')
        try
            proj = R.ProjectedCRS;
            [x, y] = projfwd(proj, lat, lon);
            [col, row] = R.worldToIntrinsic(x, y);
        catch
            [col, row] = R.worldToIntrinsic(lon, lat);
        end
    elseif isa(R, 'map.rasterref.GeographicCellsReference') || ...
           isa(R, 'map.rasterref.GeographicPostingsReference')
        [col, row] = R.geographicToIntrinsic(lat, lon);
    else
        [col, row] = R.worldToIntrinsic(lon, lat);
    end
    row = round(row); col = round(col);
end

%% ============================================================
%% NOTES: HOW TO DERIVE FLOW ACCUMULATION
%% ============================================================
%
% ---- ArcGIS Pro (recommended) ----
%   1. Fill DEM sinks:         Spatial Analyst > Hydrology > Fill
%   2. Flow direction raster:  Spatial Analyst > Hydrology > Flow Direction (D8)
%   3. Flow accumulation:      Spatial Analyst > Hydrology > Flow Accumulation
%   Save output as GeoTIFF.
%
% ---- QGIS (free) ----
%   1. Processing > SAGA > Terrain Analysis > Fill Sinks (Wang & Liu)
%   2. Processing > SAGA > Terrain Analysis > Flow Accumulation (Top-Down)
%   Save output as GeoTIFF.
%
% ---- MATLAB with TopoToolbox (free, highly recommended for geomorphology) ----
%   https://topotoolbox.wordpress.com/
%   DEM = GRIDobj('your_dem.tif');
%   DEM = fillsinks(DEM);
%   FD  = FLOWobj(DEM);
%   A   = flowacc(FD);
%   GRIDobj2geotiff(A, 'sample_flowacc.tif');
%
% ---- Threshold guidance ----
%   Flow accumulation is usually expressed as an upstream pixel count.
%   The appropriate threshold depends on raster resolution and drainage
%   geometry. Test candidate values and inspect the QA table and map;
%   lowering the threshold includes more headwater pixels.
