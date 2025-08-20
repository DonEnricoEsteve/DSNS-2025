%% Initialization
% Set path to directories
path_ft = '/Volumes/Seagate/foodEx_data/fieldtrip-20250523';
path_files = '/Volumes/Seagate/foodEx_data'; 
path_other_files = '/Volumes/Seagate/foodEx_data/DON';


paths = {path_ft, path_files, path_other_files};

missing = {};
% make sure all paths given are existing folder 
for i=1:length(paths)
    if ~isfolder(paths{i})
        missing{end+1} = paths{i};
    end
end

if ~isempty(missing)
        error('The following paths were not found:  \n %s  \n ', strjoin(missing, newline))
end

% Start fieldtrip
addpath(path_ft);
ft_defaults;

% Define subjects (according to subjects' folder names) and trigger values
subjects = {};
ranges = [3:16, 18:19, 21:22, 25:29, 31:34, 36:37, 39:41, 43:45, 47, 49:54];
for i = ranges
    subjects{end+1} = sprintf('sub_%03d', i);
end

trigVal=[8 10 12 14 20 22 24 110 112 114 120 122 124 210 212 214 220 222 224]; % corresponds to each of the 19 conditions

% Atlas for reference
atlas = ft_read_atlas(fullfile(path_ft, 'template/atlas/aal/ROI_MNI_V4.nii'));
atlas = ft_convert_units(atlas,'cm');

% Template mri
mri = ft_read_mri(fullfile(path_ft, 'external/spm8/templates/T1.nii'));

% Template grid (own)
datafile = fullfile(path_other_files, 'don_template_grid.mat');
load(datafile); % loads a variable named "template_grid"

% Neighbors
datafile = fullfile(path_other_files, 'don_neighbors_new.mat');
load(datafile); % loads a variable named "neighbors"

% Regions
datafile = fullfile(path_other_files, 'don_regions.mat');
load(datafile); % loads a variable named "regions"

%% -------------------------------------- FREQUENCY ANALYSIS SECTION --------------------------------------

%% Subtract ERF from data
function calculate_induced(path_files, subjects, trigVal)
    % Subtracts the average evoked response (ERF) from each trial of each condition, keeping only the induced response
    
    % Parameters:
    % path_files (str): directory containing all subject folders
    % subjects (cell): cell containing strings of subject folders
    % trigVal (array): Trigger values for each condition
    
    arguments
        path_files (1,1) string {mustBeFolder, mustBeNonempty}
        subjects (1, :) cell {mustBeText, mustBeNonempty}
        trigVal (1, :) double {mustBeNumeric, mustBeNonempty}
    end

    for sub = 1:length(subjects)

        fprintf("Subject: %s processing, calculate_induced  \n", subjects{sub})

        % --- Load within-subject preprocessed data ---
        datapath = fullfile(path_files, subjects{sub});
        
        % checking for path correctedness
        if ~isfolder(datapath)
            warning('The following folder does not exist: %s', datapath)
            fprintf('Ending current subject processing prior to completion  \n')
            continue   % skip to next subject in loop
        end

        datafile = fullfile(datapath, "ICA_dataorig.mat");
        if ~isfile(datafile)
            warning('ICA_dataorig.mat not found')
            fprintf('Ending current subject processing prior to completion  \n')
            continue   % skip to next subject in loop
        end
        
        % Load
        load(datafile, 'ICA_dataorig') % loads a variable named "ICA_dataorig (bpfreq = [1 90]) (1x1 struct)"
        
        % Verify correct loaded variable structure
        if ~isstruct(ICA_dataorig)
            warning("Variable ICA_dataorig must be a structure with fields:  " + ...
                "time - epoched time of trials (1 x trials cell (1 x no. samples double))  "+ ...
                "trial - data for each trial (1 x trials cell (no. channels x no. samples double))  "+ ...
                "label - channel names (1 x no. channels double)  " + ...
                "grad -   " + ...
                "trialinfo - trial trigger values (conditions) in the 1st column (double) " + ...
                "cfg -  ")
            fprintf('Ending current subject processing prior to completion  \n')
            continue

        elseif ~isfield(ICA_dataorig, {'trialinfo', ' trial', 'time'})
            warning("ICA_dataorig is missing one or all the following fields:" + ...
                "time - epoched time of trials (1 x trials cell (1 x no. samples double))  "+ ...
                "trial - data for each trial (1 x trials cell (no. channels x no. samples double))  " + ...
                "trialinfo - trial trigger values (conditions) in the 1st column (double) ")         
            fprintf('Ending current subject processing prior to completion  ')
            continue

        elseif ~isfield(ICA_dataorig, 'grad') || ~isfield(ICA_dataorig.grad, 'unit')
            warning("ICA_dataorig must have grad.unit field")
            fprintf('Ending current subject processing prior to completion \n ') 
            continue
        
        elseif ~iscell(ICA_dataorig.time) || ~iscell(ICA_dataorig.trial)
            warning("The following fields in ICA_dataorig should be of type cell:" + ...
            "time - epoched time of trials (1 x trials cell (1 x no. samples double))  " + ...
            "trial - data for each trial (1 x trials cell (no. channels x no. samples double))  ")
            fprintf('Ending current subject processing prior to completion \n ') 
            continue

        elseif ~(length(ICA_dataorig.trialinfo(:,1)) == length(ICA_dataorig.trial))
            warning("ICA_dataorig.trialinfo(:,1) indicating the trigger value of each trial  " + ...
                "and ICA_dataorig.trial containing the data for each trial must have the same length  " + ...
                "(1 x no. trials)  ")
            fprintf('Ending current subject processing prior to completion \n ')
            continue

        elseif ~all(ismember(unique(ICA_dataorig.trialinfo(:,1)), unique(trigVal)))
            warning("trigVal should have the unique trigger values as specified for the trials in  " + ...
                "ICA_dataorig.trialinfo(:,1)  ")
            fprintf('Ending current subject processing prior to completion \n ')
            continue
        end
        
        % Check that number of time points recorded matches number of data
        % samples
        nTime = cellfun(@length, ICA_dataorig.time);
        nSample = cellfun(@(x) size(x,2), ICA_dataorig.trial);
        
        mismatched = find(nTime ~= nSample);
            
        if ~isempty(mismatched)
            warning("Sample mismatch in between ICA_dataorig.time and ICA_dataorig.trial " + ...
                "in the following trials: %s", strjoin(string(mismatched), ', '));
        end

        % --- Load within-subject ERF ---
        datafile = fullfile(datapath, "ERF_within", "don_ERF.mat");
        
        % Checking for path correctedness
        if ~isfile(datafile)
            warning('%s not found', datafile);
            fprintf('Ending current subject processing prior to completion \n ')
            continue   % skip to next subject
        end
        
        % Load
        load(datafile, 'ERF') % loads a variable named "ERF" (1x19 cell)
        
        % Verify correct loaded variable structure
        if ~iscell(ERF)
            warning("ERF variable must be a cell of size 1 x no. conditions, each cell  " + ...
                "containing the following fields:  " + ...
                " ")
            fprintf('Ending current subject processing prior to completion \n ')
            continue
        elseif ~(length(ERF)<=length(trigVal))
            warning("ERF can't have more conditions than listed in trigVal, " + ...
                "trigVal must list all possible conditions in order")
            fprintf('Ending current subject processing prior to completion \n ')
            continue
        end

        skipOuter = false;

        % Verify correct loaded variable, internal cell structure
        for  cond=1:length(ERF)
            if isstruct(ERF{1, cond}) 
                if ~isfield(ERF{1, cond}, 'avg')
                     warning("ERF should contain cells as the number of conditions " + ...
                    "with each cell being a structrue with field 'avg'.   Cell %d doesn't have fiels 'avg'.", ...
                    con)
                     fprintf('Ending current subject processing prior to completion \n ')
                     skipOuter = true;
                     break; 
                end
            else 
                warning("ERF should contain cells as the number of conditions " + ...
                    "with each cell being a structrue with field 'avg'.   Cell %d is not a structure.", ...
                    con)
                fprintf('Ending current subject processing prior to completion \n ')
                skipOuter = true;
                break;
            end
        end
        
        if skipOuter
            continue
        end

        % --- Initialize a vector to contain preprocessed data without ERF ---
        dataorig_minus_ERF = ICA_dataorig;
        

        % --- Calculate within-subject induced data ---
        
        skipped_con = 0;
        % For each condition:
        for con = 1:length(trigVal)
    
            % Take the indices of each condition relative to the trial
            trl = find(ICA_dataorig.trialinfo(:,1) == trigVal(con));
            
            % Verify trl is not empty - there were trials found with the
            % specific trigger value
            if isempty(trl)
                warning('No trials of condition %d were found.', trigVal(con))
                fprintf('Ending current subject processing prior to completion \n ')
                skipped_con = skipped_con + 1;
                continue
            end

            % Take the values within the indices from the .trl field
            con_trl_vals = ICA_dataorig.trial(:, trl);
            
            % Subtract the avg field of ERF from each subsetted trial of ICA_dataorig
            for k = 1:numel(con_trl_vals)
                trial = con_trl_vals{k};
                cond_ERF = con - skipped_con;
                evoked = ERF{1, cond_ERF}.avg;
                if ~isequal(size(trial), size(evoked))
                    warning("Condition: %d. Trial number %d and ERF cell %d, don't have same dimensions.", trigVal(con), k,  cond_ERF)
                    fprintf("Moving to the next trial of the same condition. \n ")
                    continue
                end
                dataorig_minus_ERF.trial{trl(k)} = trial - evoked;
            end
        end
    
        % Save preprocessed data without ERF
        savefile = fullfile(datapath, 'dataorig_minus_ERF.mat');
        save(savefile, 'dataorig_minus_ERF'); % yields a 1x1 struct
        
        fprintf("Subject %s was processed successfuly \n ", subjects{sub})
    end
end

calculate_induced(path_files, subjects, trigVal)

%% Frequency Analysis
function calculate_freq(path_files, subjects, trigVal, time_interval, baseline)
    % Parameters:
    % path_files (str): Directory containing all subject folders
    % subjects (cell): Cell containing strings of subject folders
    % trigVal (array): Trigger values for each condition
    % time_interval (array): Range of data truncation
    % baseline (bool): is time_interval baseline, default = false
    
    arguments
        path_files (1,1) string {mustBeFolder, mustBeNonempty}
        subjects (1, :) cell {mustBeText, mustBeNonempty}
        trigVal (1, :) double {mustBeNumeric, mustBeNonempty}
        time_interval (1, 2) double {mustBeNumeric, mustBeNonempty}
        baseline (1,1) logical = false
    end
    
    for sub = 1:length(subjects)
                
        fprintf("Subject: %s processing, calculate_freq  \n", subjects{sub})        
        
        % --- Load within-subject preprocessed data ---
        datapath = fullfile(path_files, subjects{sub});
        datafile = fullfile(datapath, "dataorig_minus_ERF.mat");

        % Verify path exists
        if ~isfile(datafile)
            warning("%s does not exist", datafile)
            fprintf('Ending current subject processing prior to completion \n ')
            continue
        end
        
        % Load
        load(datafile, 'dataorig_minus_ERF') % loads a variable named "dataorig_minus_ERF"
        
        % Convert grad units
        dataorig_minus_ERF.grad = ft_convert_units(dataorig_minus_ERF.grad, 'cm');
        
        leave_idx = find(dataorig_minus_ERF.time{1}<=max(time_interval) & dataorig_minus_ERF.time{1}>=min(time_interval));
        
        if isempty(leave_idx)
            warning("The desired time_interval was not found in dataorig_minus_ERF.time cells")
            fprintf('Ending current subject processing prior to completion \n ')
        end
        
        % --- Truncate data and time vector to time_interval range ---

        % For all trials
        for i = 1:length(dataorig_minus_ERF.trial)

            % Leave only time_interval data
            dataorig_minus_ERF.trial{i} = dataorig_minus_ERF.trial{i}(:, leave_idx); 

            % Apply the same truncation to the corresponding time vector
            dataorig_minus_ERF.time{i} = dataorig_minus_ERF.time{i}(leave_idx);
        end
         
        
        % Check if desired freq resolution is better than maximal
        skipped_con = 0;
        freq_resul = 2;
        fs = 1017;
        n_samples = length(dataorig_minus_ERF.time{1});
        freq_resul_max = fs/n_samples;

        if freq_resul_max>freq_resul
            warning("Desired frequency resolution for freq analysis is better than the  " + ...
                "inherent signal resolution due to low n_samples, ft functions will interpolate the non-existing  " + ...
                "frequencies.")
        end
        
         % --- Calculate frequency representation per condition onn truncated data ---

        for c = 1:length(trigVal)
            cfg = [];

            trial_idx = find(dataorig_minus_ERF.trialinfo(:,1) == trigVal(c)); 

            if isempty(trial_idx)
                skipped_con = skipped_con +1;
                warning("Trials that correspond to condition %d were not found.   Skipping to next condition " + ...
                    "for frequency representation calculation", trigVal(c))
                continue
            end

            cfg.trials = find(dataorig_minus_ERF.trialinfo(:,1) == trigVal(c));
            cfg.output = 'pow';   % Compute CSD matrices (powandcsd for channel connectivity)
            cfg.pad = 'nextpow2'; % for more efficient FFT 

            % The padding will determine your spectral resolution. 
            % cfg.keeptrials = 'yes'; % phase estimates are required for each individual trial (required for fourier)
            
            % For frequency analysis
            cfg.method = 'mtmfft';  % assumes stationarity
            cfg.taper = 'hanning';  % no smoothing, if required, use 'dpss' and enable tapsmofrq
            cfg.foi = freq_resul:freq_resul:90;       % frequency resolution is determined by time window length (1/800 ms=1.25 Hz Nyquist)
            % cfg.tapsmofrq = 2;

            freq{1, c-skipped_con} = ft_freqanalysis(cfg, dataorig_minus_ERF);
        end
    
        % --- Save ---

        % Check path for saving exists, if not, make directory
        if ~isfolder(fullfile(datapath, 'TFR_within'))
            warning("The path %s for saving the frequency data, does not exist. Creating the folder.", fullfile(datapath, 'TFR_within'))
            mkdir(fullfile(datapath, 'TFR_within'))
        end

        if baseline == true
            savefile = fullfile(datapath, 'TFR_within', 'don_freq_allconsbaseline.mat');
        else
            savefile = fullfile(datapath, 'TFR_within', 'don_freq_allcons.mat');
        end
        save(savefile, 'freq'); % yields a 1x1 struct
        fprintf("Subject %s was processed successfuly \n ", subjects{sub})
    end
    
end
    
tic
time_interval = [0.0, 0.8];
calculate_freq(path_files, subjects, trigVal, time_interval, false)
toc 
% Elapsed time is 1605.629980 seconds.

%% Concatenate freqs across subjects
function contatenate_freqs(path_files, subjects, load_file, save_file)
    
    % Initialize cell
    freq_across = {};
   
    % --- Concatenate powerspectrum per consition across subjects ---

    % Load within-subject ERF
    for sub = 1:length(subjects)
        fprintf("Subject: %s processing, contatenate_freqs  \n", subjects{sub})  
        datapath = fullfile(path_files, subjects{sub});
        datafile = fullfile(datapath, "TFR_within", load_file);

        % Verify path exists
        if ~isfile(datafile)
            warning("%s does not exist", datafile)
            fprintf('Ending current subject processing prior to completion \n ')
            continue
        end

        load(datafile, 'freq') % loads a variable named "freq"

        % For each condition  
        for con = 1:size(freq,2)
            % Contatenate the freqs without RMS. Used for plotting.
            freq_across{sub, con} = freq{1, con};
        end
        fprintf("Subject %s was concatenated successfuly \n ", subjects{sub})
    end

    % Save the contatenated cell 
    datapath = fullfile(path_files, "TFR_DSNS");
    
    % Verify path exists
    if ~isfolder(datapath)
        error("%s does not exist", datapath)
    end
    savefile = fullfile(datapath, save_file);
    save(savefile, 'freq_across', '-v7.3'); % yields a 42x2 cell each containing a 1x1 struct
    fprintf("All subjects were concatenated successfully")
end

tic
contatenate_freqs(path_files, subjects, "don_freq_allcons.mat", "don_freq_across_allcons.mat")
contatenate_freqs(path_files, subjects, "don_freq_allconsbaseline.mat", "don_freq_across_allconsbaseline.mat")
toc 
% Elapsed time is 176.489619 seconds.

%% Load concatenated freq data for statistics and plotting
datapath = fullfile(path_files, "TFR_DSNS");
datafile = fullfile(datapath, "don_freq_across_allcons.mat");
freq_across = load(datafile, 'freq_across').freq_across; % loads a variable named "freq_across"

datapath = fullfile(path_files, "TFR_DSNS");
datafile = fullfile(datapath, "don_freq_across_allconsbaseline.mat");
freqbase_across = load(datafile, 'freq_across').freq_across; % loads a variable named "freq_across"

%% Normalize and extract powspctrm values from frequency data
% TO BE USED FOR REPEATED MEASURES ANOVA IN JAMOVI
function avgPower = extract_2D_power(baseline_cell, poststim_cell, channels, freqrange)
    
    arguments
        baseline_cell cell {mustBeNonempty}
        poststim_cell cell {mustBeNonempty}
        channels {mustBeText, mustBeNonempty} 
        freqrange (1,2) double {mustBeNumeric, mustBeNonempty}
    end

    all_cells_size = cellfun(@(x) size(x.powspctrm), [baseline_cell(:); poststim_cell(:)], 'UniformOutput', false);

    if ~isequal(size(baseline_cell), size(poststim_cell))
        error("Length of baseline and post-stimulus cells should be equal")
    elseif ~all(all(cellfun(@(x) isfield(x, 'label'), baseline_cell))) ||  ~all(all(cellfun(@(x) isfield(x, 'label'), poststim_cell)))...
        ||  ~all(all(cellfun(@(x) isfield(x, 'freq'), baseline_cell))) ||  ~all(all(cellfun(@(x) isfield(x, 'freq'), poststim_cell)))...
        || ~all(all(cellfun(@(x) isfield(x, 'powspctrm'), baseline_cell))) ||  ~all(all(cellfun(@(x) isfield(x, 'powspctrm'), poststim_cell)))
        error("poststim_cell or baseline_cell lack the following fields in one or more cells: label, freq, powspctrm.")
    elseif ~(channels == "all" || (all(cellfun(@(x) ismember(channels, x.label), baseline_cell)) && all(cellfun(@(x) ismember(channels, x.label), poststim_cell))))
        error('channels should be "all" or consist of a list of channels found in baseline_cell.label and poststim_cell.label.')
    elseif ~(all(all(cellfun(@(x) isequal(length(x.freq),size(x.powspctrm,2)), baseline_cell))) && all(all(cellfun(@(x) isequal(length(x.freq),size(x.powspctrm,2)), poststim_cell))))
        error("All cells in baseline_cell and poststim_cell should have an equal number of frequencies in  freq field and columns of powspctrm field.")
    elseif ~(all(all(cellfun(@(x) isequal(length(x.label),size(x.powspctrm,1)), baseline_cell))) && all(all(cellfun(@(x) isequal(length(x.label),size(x.powspctrm,1)), poststim_cell))))
        error("All cells in baseline_cell and poststim_cell should have an equal number of channels in label field and rows of powspctrm field.")
    elseif ~isequal(all_cells_size{:})
        error("All cells should have an equal powspctrm, in baseline_cell and poststim_cell")
    end

    [nSubj, nCond] = size(poststim_cell);
    avgPower = nan(nSubj, nCond); % initialize avgPower matrix with NaNs

    for s = 1:nSubj
        for c = 1:nCond
            baseData = baseline_cell{s, c};
            postData = poststim_cell{s, c};

            % --- Select desired channels ---
            cfg = [];
            cfg.channel = channels;
            baseSel = ft_selectdata(cfg, baseData);
            postSel = ft_selectdata(cfg, postData);

            % --- Normalize per channel × frequency (dB) ---
            normPow = 10 * log10(postSel.powspctrm ./ baseSel.powspctrm); 

            % --- Average across channels ---
            chanAvg = mean(normPow, 1, 'omitnan'); % 1 × nFreq

            % --- Average within the frequency band ---
            fidx = postSel.freq >= freqrange(1) & postSel.freq <= freqrange(2);
            if isempty(fidx)
                error("Could find no data for the desired freqrange in cell {%d, %d}.", s ,c)
            end
            avgPower(s, c) = mean(chanAvg(fidx), 'omitnan');
        end
    end
end
tic
avg2DPower = extract_2D_power(freqbase_across, freq_across, "all", [25 50]);
toc

if ~isequal(size(freq_across), size(avg2DPower))
    error('avg2DPoer should be equal in size to freq_across/freqbase_across (no_subjects x no_conditions)')
end

% Global (all channels)
% θ-band (4–8 Hz) 
% α-band (8–13 Hz) 
% β-band (13–25 Hz)
% γ-band (25–50 Hz)

%% Plot across subject frequency spectrum across channels
function plot_freq_across_groups(freq_across, freq_baseline, condGroups, channels, tle)
% freq_across:  nSubj × nCon cell of freq structs
% freq_baseline: same size as freq_across
% condGroups:    cell array of vectors, e.g., {[2:4], [5:7], [8:10]}
% channels:      channel selection for averaging (e.g., 'all' or {'MEG'})
% tle:           plot title
    
    arguments
        freq_across cell {mustBeNonempty}
        freq_baseline cell {mustBeNonempty}
        condGroups cell {mustBeNonempty}
        channels {mustBeText, mustBeNonempty} 
        tle {mustBeTextScalar}
    end
    
    % Verify valid input for condGroups
    
    conditions = [];
    for group=condGroups
        conditions = [conditions group{1}];
    end
        
    if ~all(cellfun(@(x) isnumeric(x) && ~isempty(x) && all(mod(x,1)==0), condGroups))
        error("All cells in condGroups must contain integer arrays")
    elseif ~all(ismember(conditions, 1:size(freq_across, 2)))
        error("The conditions in condGroups are out of freq_across and freq_baseline data condition index range.")
    elseif ~isequal(conditions, unique(conditions))
        error("There are repeating conditions in condGroups")
    end

    nSubj = size(freq_across, 1);
    nGroups = numel(condGroups);

    % Preallocate
    avg_freq = cell(nGroups, 1);
    for g = 1:nGroups
        avg_freq{g} = cell(nSubj, 1);
    end

    % --- Loop over subjects ---
    for s = 1:nSubj
        for g = 1:nGroups
            conds = condGroups{g};
            fprintf("Averaging across conditions %d, subject %d - plot_freq_across_groups", nSubj, conds)
            normConds = cell(numel(conds), 1);

            % --- Loop over conditions in this group ---
            for k = 1:numel(conds)
                postData = freq_across{s, conds(k)};
                baseData = freq_baseline{s, conds(k)};

                % Select channels
                cfg = [];
                cfg.channel = channels;
                postSel = ft_selectdata(cfg, postData);
                baseSel = ft_selectdata(cfg, baseData);

                % dB normalization per channel × frequency
                if ~(all(all(postSel.powspctrm>0)) &&  all(all(baseSel.powspctrm>0)))
                    error("All freq_across and freq_baseline powspectrm should have positive values. \n" + ...
                        "The following cell doesn't have positive values: {%d, %d}", s, conds(k))
                end

                normPow = 10 * log10(postSel.powspctrm ./ baseSel.powspctrm);

                % Average across channels
                chanMean = mean(normPow, 1, 'omitnan');
                
                if ~isequal(size(chanMean), size(postSel.freq))
                    error("chanMean size is (%d, %d) and doesn't equal to the number of frequencies in postSel.freq (%d)", size(chanMean,1), size(chanMean,2), length(postSel.freq))
                end

                % Create freq struct for this normalized condition
                normStruct = postSel;
                normStruct.powspctrm = chanMean;
                normStruct.label = {'avg'};
                normStruct.dimord = 'chan_freq';
                normConds{k} = normStruct;
            end

            % Average across all conditions in this group for this subject
            cfg = [];
            cfg.keepindividual = 'no';
            cfg.parameter = 'powspctrm';
            avg_freq{g}{s} = ft_freqgrandaverage(cfg, normConds{:});
        end
    end
    fprintf("Processed all subjects and all condGroups")

    % --- Grand average across subjects for each group ---
    ga_freq = cell(nGroups, 1);
    cfg = [];
    cfg.parameter = 'powspctrm';
    for g = 1:nGroups
        ga_freq{g} = ft_freqgrandaverage(cfg, avg_freq{g}{:});
    end

    % --- Plot ---
    figure('Position', [100 100 600 400]); 
    hold on;

    % Define some colors and line styles for multiple groups
    colors = [0 0 0; 0 0 1; 0 0.6 0; 0.9 0.6 0; 0.6 0 0.6; 0.7 0.7 0.7];
    styles = {'-', '--', ':', '-.', '-', '--'}; % will cycle if >4 groups

    for g = 1:nGroups
        plot(ga_freq{g}.freq, squeeze(ga_freq{g}.powspctrm), ...
             'LineWidth', 2, ...
             'Color', colors(g,:), ...
             'LineStyle', styles{mod(g-1,numel(styles))+1});
    end

    % --- Mark effects ---
    xlim([0 90]);
    ylim([-5 1]);
    y_bounds = ylim; 
    xranges = [8 13; 25 50]; % e.g., beta band

    for i = 1:size(xranges, 1)
        fill([xranges(i,1) xranges(i,2) xranges(i,2) xranges(i,1)], ...
             [y_bounds(1) y_bounds(1) y_bounds(2) y_bounds(2)], ...
             [0.9 0.9 0.9], ...
             'EdgeColor', 'k', ...
             'LineWidth', 0.5, ...
             'LineStyle', '--', ...
             'FaceAlpha', 0.5);
    end

    yline(0, 'r', 'LineWidth', 2);

    % Create dynamic legend
    customNames = {'Food-Novel', 'Food-Repeated', 'Positive-Novel', 'Positive-Repeated', ...
        'Neutral-Novel', 'Neutral-Repeated'};
    
    if nGroups > numel(customNames)
        warning("The number of condGroups for analysis is larger than the number of customNames for the plot legend. \n CustomNames might not match the condition groups.")
    elseif nGroups < numel(customNames)
        warning("The number of condition groups is smaller than the number of customNames, customNames might not match the condition groups.")
    end

    nGroups = min(nGroups, numel(customNames)); % if the number of groups is larger than legend handles, 
    % use the length of legend handles and if otherwise use the number of condGroups 

    lgdNames = arrayfun(@(g) sprintf('%s', customNames{g}), 1:nGroups, 'UniformOutput', false);
    lgd = legend(lgdNames);
    title(lgd, 'Category-Presentation'); 

    xlabel('Frequency (Hz)', 'FontSize', 14);
    ylabel('Global Power Change (dB)', 'FontSize', 14);
    title(tle, 'FontSize', 20);
    set(gca, 'FontSize', 12);   
    grid off;
end

tic
condGroups = {2:4, 5:7, 8:10, 11:13, 14:16, 17:19};
plot_freq_across_groups(freq_across, freqbase_across, condGroups, "all", ...
    "PSD across channels and subjects: Comparison among Category-Presentation");
toc

%% -------------------------------------- TIME-FREQUENCY ANALYSIS SECTION --------------------------------------

%% Time-Frequency Analysis
function calculate_TFR(path_files, subjects, trigVal)
    
    arguments
        path_files (1,1) string {mustBeFolder, mustBeNonempty}
        subjects (1, :) cell {mustBeText, mustBeNonempty}
        trigVal (1, :) double {mustBeNumeric, mustBeNonempty}
    end

    for sub = 1:length(subjects)
        fprintf("Processing subject %d", sub)
        % Load within-subject preprocessed data
        datapath = fullfile(path_files, subjects{sub});
        datafile = fullfile(datapath, "dataorig_minus_ERF.mat");
        load(datafile) % loads a variable named "dataorig_minus_ERF"

        % Convert grad units
        dataorig_minus_ERF.grad = ft_convert_units(dataorig_minus_ERF.grad, 'cm');

        % Calculate time-frequency representation
        for c = 1:length(trigVal)
            cfg = [];
            cfg.trials = find(dataorig_minus_ERF.trialinfo(:,1) == trigVal(c));
            cfg.output = 'pow';   % Compute CSD matrices (powandcsd for channel connectivity)

            % For time-frequency analysis
            cfg.method     = 'mtmconvol';
            cfg.taper      = 'hanning';
            cfg.foi        = 4:2:90;        % start at 4 Hz
            cfg.t_ftimwin  = 3 ./ cfg.foi;  % 3 cycles
            cfg.pad        = 2;             % 2 s padding
            cfg.padtype    = 'mirror';      % safer than zero padding
            cfg.toi        = -0.3:0.01:0.8; % include baseline
            
            % error if the time range cfg.toi isn't found in data
            if ~all(cellfun(@(x) nnz(cfg.toi(1)<x & cfg.toi(end)>x), dataorig_minus_ERF.time)) 
                error("The desired time range for time-frequency analysis (cfg.toi) doesn't match the time in dataorig_minus_ERF.time.")
            end

            % NOTES
            % For foi: The freq resolution is determined by min freq and number of cycles
            % For foi: We are starting at 4 Hz since it yields a time window that safely fits the trial length of 1.1s
            % For t_ftimwin: N_cycles / f = window length (s); higher N_cycles, better frequency resolution
            % For toi: The lowest SAFE toi is (N_cycles/min_foi) = 3/4 = 0.75 s -> half window = 0.375 s
            % For pad: Pad to at least trial length + longest window length = 1.1 s + 0.75 s = 1.85 s ~ 2 s safe padding
            
            % “trustworthy” ERD/ERS range is roughly 0.075 to 0.425 s for the lowest frequency
            % For toi: toi_min = start_toi + 0.25 | toi_max = end_toi - 0.25 = [-0.3+0.25= -0.05s, 0.8-0.25= 0.55s]

            TFR{1, c} = ft_freqanalysis(cfg, dataorig_minus_ERF);
        end
    
        % Save
        savefile = fullfile(datapath, 'TFR_within', 'don_newTFR_allconsbp.mat');
        save(savefile, 'TFR'); % yields a 1x1 struct
        fprintf("Finished processing subject %s", sub)

    end
end
    
tic
calculate_TFR(path_files, subjects, trigVal)
toc 
% Elapsed time is 6234.724326 seconds.

%% Normalize and concatenate TFRs across subjects
function normalize_and_concatenate_TFRs(path_files, subjects)
    
    arguments
        path_files (1,1) string {mustBeFolder, mustBeNonempty}
        subjects (1, :) cell {mustBeText, mustBeNonempty}
    end

    nSubj = length(subjects);
    normTFR_across = {}; % will be nSubj x nConditions

    for sub = 1:nSubj
        fprintf("Processing subject %d, normalize_and_concatenate_TFRs", sub)
        % --- Load within-subject preprocessed data ---
        datapath = fullfile(path_files, subjects{sub});
        datafile = fullfile(datapath, "TFR_within", "don_newTFR_allconsbp.mat");
        load(datafile, 'TFR')  % loads cell array of TFRs

        % Preallocate per-subject normalized TFR
        nCon = size(TFR,2);
        normTFR = cell(1,nCon);

        % --- Baseline normalization per condition ---
        for c = 1:nCon
            cfg = [];
            cfg.baseline    = [-0.3 0]; % baseline interval
            cfg.baselinetype = 'db';    % dB: 10*log10(post/baseline)
            cfg.parameter    = 'powspctrm';
            normTFR{1,c} = ft_freqbaseline(cfg, TFR{1,c});

            % Store normalized TFR into across-subject cell
            normTFR_across{sub, c} = normTFR{1,c};
        end
    end

    % --- Save the concatenated cell across subjects ---
    datafile = fullfile(path_files, "TFR_DSNS");
    savefile = fullfile(datafile, "don_normTFR_across_allconsbp.mat");
    save(savefile, 'normTFR_across', '-v7.3'); % nSubj x nCon cell array
    fprintf("Finished processing all subjects and saving")
end

normalize_and_concatenate_TFRs(path_files, subjects)

%% Load concatenated TFR data for statistics and plotting
datapath = fullfile(path_files, "TFR_DSNS");
datafile = fullfile(datapath, "don_normTFR_across_allconsbp.mat");
load(datafile) % loads a variable named "normTFR_across"

%% Extract powspctrm values from time-frequency data
% TO BE USED FOR REPEATED MEASURES ANOVA IN JAMOVI
function avgPower = extract_3D_power(data_across, channels, freqrange, timerange)

    arguments
        data_across cell 
        channels  {mustBeText, mustBeNonempty} 
        freqrange (1,2) double {mustBeNumeric, mustBeNonempty}
        timerange (1,2) double {mustBeNumeric, mustBeNonempty}
    end

    if ~all(cellfun(@(x) nnz(timerange(1)<x.time & timerange(end)>x.time), data_across))
         error("The desired time range doesn't match the time in data_across.")

    elseif ~all(cellfun(@(x) nnz(freqrange(1)<x.freq & freqrange(end)>x.freq), data_across)) 
         error("The desired frequency range doesn't match the frequencies in data_across.")

    elseif ~all(cellfun(@(x) isequal(size(x.powspctrm), [length(x.label), length(x.freq), length(x.time)]), data_across))
    error("data_across.powspctrm should have the size of (channels, frequencies, time).")
    end

    [nSubj, nCond] = size(data_across);
    avgPower = nan(nSubj, nCond);  % Preallocatear

    for s = 1:nSubj
        fprintf("Processing subject %d, extract_3D_power \n", s)
        for c = 1:nCond
            % Extract single condition struct
            freqData = data_across{s, c};
            % Step 1: Select desired channels
            cfg = [];
            cfg.channel = channels;
            selData = ft_selectdata(cfg, freqData);

            % Average over channels
            chanAvg = squeeze(mean(selData.powspctrm, 1, 'omitnan'));  % freq × time

            % Step 2: Select desired frequency range
            fidx = selData.freq >= freqrange(1) & selData.freq <= freqrange(2);
            freqAvg = mean(chanAvg(fidx, :), 1, 'omitnan');  % 1 × time

            % Step 3: Select desired time range
            tidx = selData.time >= timerange(1) & selData.time <= timerange(2);
            timeAvg = mean(freqAvg(tidx), 'omitnan');  % scalar

            if length(timeAvg)>1 || isempty(timeAvg)
                error("timeAvg should be a scalar, average time-frequency data for subject and condition.")
            end

            avgPower(s, c) = timeAvg;
        end
    end
end

tic
avg3DPower = extract_3D_power(normTFR_across, "all", [4 8], [0.100 0.200]);
toc

%% Calculate frequency sensor-level cluster statistics
function TFR_cluster_statistics(normTFR_across, path_files, neighbors, freqrange, timerange, con1, con2)

    % Average selected conditions within each subject
    nSubj = size(normTFR_across, 1);
    data1 = cell(nSubj, 1); % con1
    data2 = cell(nSubj, 1); % con2

    for s = 1:nSubj
        cfg = [];
        cfg.keepindividual = 'no'; % average within subject
        cfg.parameter = 'powspctrm';

        data1{s} = ft_freqgrandaverage(cfg, normTFR_across{s, con1});
        data2{s} = ft_freqgrandaverage(cfg, normTFR_across{s, con2});
    end

    % ----- Cluster-based statistics configuration -----
    cfg = [];
    cfg.frequency     = freqrange;        % e.g., [8 12] or 'all'
    cfg.latency       = timerange;        % e.g., [0.2 0.6] or 'all'
    cfg.channel       = {'MEG', '-A17', '-A203'};
    cfg.parameter     = 'powspctrm';
    cfg.method        = 'montecarlo';
    cfg.statistic     = 'depsamplesT';    % within-subject design
    cfg.correctm      = 'cluster';
    cfg.neighbours    = neighbors;

    cfg.clusteralpha      = 0.05;
    cfg.tail              = 0;
    cfg.clustertail       = 0;
    cfg.alpha             = 0.025;
    cfg.numrandomization  = 1000;

    cfg.minnbchan         = 2;        % Ignore clusters with <2 neighbors
    cfg.correcttail       = 'prob';   % Optional: two-sided cluster correction

    % ----- Design matrix -----
    nsub = nSubj;
    design = zeros(2, 2*nsub);
    design(1,:) = [1:nsub 1:nsub];                % subject indices
    design(2,:) = [ones(1,nsub) 2*ones(1,nsub)];  % condition labels

    cfg.design = design;
    cfg.uvar   = 1;   % subject
    cfg.ivar   = 2;   % condition

    % ----- Compute cluster statistics (Con2 vs Con1) -----
    % Swap order: first all con2, then all con1
    TFR_stat = ft_freqstatistics(cfg, data2{:}, data1{:});

    % ----- Save result -----
    datapath = fullfile(path_files, "TFR_DSNS/FreqStats_New");
    savefile = fullfile(datapath, "don_normTFR_stat_allfreq_alltime_FS2vFS1.mat");
    save(savefile, 'TFR_stat', '-v7.3');

end

tic
TFR_cluster_statistics(normTFR_across, path_files, neighbors, "all", "all", ...
    2, 5)
toc

%Elapsed time is 3645.774012 seconds.

%% Load TFR stats data for plotting
datapath = fullfile(path_files, "TFR_DSNS/FreqStats_New");
datafile = fullfile(datapath, "don_normTFR_stat_allfreq_alltime_2v1.mat");
load(datafile) % loads a variable named "TFR_stat"

%% Plot TFR statistics clusters
function plot_TFR_cluster_statistics(TFR_stat, time, freq, avgfreq)
    arguments
        TFR_stat
        time
        freq
        avgfreq
    end
    cfg = [];
    cfg.latency = time;
    cfg.frequency = freq;

    if ~isempty(avgfreq)
        cfg.avgoverfreq = 'yes';
    else
        cfg.avgovertime = 'yes';
    end
      
    TFR_stat_selected = ft_selectdata(cfg, TFR_stat);
    
    cfg = [];
    cfg.layout = '4D248_helmet.mat';
    cfg.alpha = 0.05;
    cfg.subplotsize = [3 3];
    % cfg.colorbar = 'yes';
    cfg.zlim = [-3 3];
    cfg.parameter= 'stat';
    ft_clusterplot(cfg, TFR_stat_selected);
end

% Keep last argument empty to avg. over time, otherwise avg. over freq
plot_TFR_cluster_statistics(TFR_stat, [0.3 0.5], [8 25], 'avgfreq')

%% Plot across-subject TFR spectrum across channels with 3 subplots
function plot_TFR_across_single(normTFR_across, con1, con2, channels, tle, figtitle)
arguments
    normTFR_across
    con1
    con2
    channels
    tle
    figtitle
end
    % tle: a 1x3 cell array for the subplot titles {Con1, Con2, Diff}

    nSubj = size(normTFR_across, 1);
    avg_TFR1 = cell(nSubj, 1);
    avg_TFR2 = cell(nSubj, 1);

    % ----- Loop through subjects -----
    for s = 1:nSubj
        % Average over selected conditions within subject
        cfg = [];
        cfg.keepindividual = 'no';
        cfg.parameter = 'powspctrm';
        data1 = ft_freqgrandaverage(cfg, normTFR_across{s, con1});
        data2 = ft_freqgrandaverage(cfg, normTFR_across{s, con2});

        % Select channels and average across them
        cfg = [];
        cfg.channel = channels;
        tmp1 = ft_selectdata(cfg, data1);
        tmp2 = ft_selectdata(cfg, data2);

        % Average over channels
        tmp1.powspctrm = mean(tmp1.powspctrm, 1); % freq x time
        tmp2.powspctrm = mean(tmp2.powspctrm, 1);

        % Update labels and dimord for 1 channel
        tmp1.label = {'avg'}; 
        tmp2.label = {'avg'};
        tmp1.dimord = 'chan_freq_time';
        tmp2.dimord = 'chan_freq_time';

        avg_TFR1{s} = tmp1;
        avg_TFR2{s} = tmp2;
    end

    % ----- Grand average across subjects -----
    cfg = [];
    cfg.parameter = 'powspctrm';
    ga_TFR1 = ft_freqgrandaverage(cfg, avg_TFR1{:});
    ga_TFR2 = ft_freqgrandaverage(cfg, avg_TFR2{:});

    % Compute difference (Con2 - Con1)
    cfg = [];
    cfg.parameter = 'powspctrm';
    cfg.operation = 'x2 - x1';
    TFR_diff = ft_math(cfg, ga_TFR1, ga_TFR2);

    % ----- Prepare figure -----
    figure('Name',figtitle,'Position',[100 100 1800 500]);

    % Data arrays for plotting
    TFRs = {ga_TFR1, ga_TFR2, TFR_diff};
    climVals = [-0.08 0.08]; % dB range

    for i = 1:3
        subplot(1,3,i);
        imagesc(TFRs{i}.time, TFRs{i}.freq, squeeze(TFRs{i}.powspctrm));
        axis xy;
        xlabel('Time (s)', 'FontSize', 14);
        ylabel('Frequency (Hz)', 'FontSize', 14);
        title(tle{i}, 'FontSize', 14);
        set(gca,'FontSize',12);
        colorbar;
        clim(climVals);

        % Optional: Draw example rectangles
        hold on;
        x_ranges = [0.3 0.5; 0.3 0.5];
        y_ranges = [8 13; 13 25];      
        for r = 1:size(x_ranges,1)
            x = [x_ranges(r,1), x_ranges(r,2), x_ranges(r,2), x_ranges(r,1), x_ranges(r,1)];
            y = [y_ranges(r,1), y_ranges(r,1), y_ranges(r,2), y_ranges(r,2), y_ranges(r,1)];
            plot(x, y, 'k--', 'LineWidth', 1);
        end
        hold off;
    end

    % Adjust colorbar labels
    subplot(1,3,1); cb1 = colorbar; ylabel(cb1, 'Global Power Change (dB)', 'FontSize', 12);
    subplot(1,3,2); cb2 = colorbar; ylabel(cb2, 'Global Power Change (dB)', 'FontSize', 12);
    subplot(1,3,3); cb3 = colorbar; ylabel(cb3, 'Difference in Global Power Change (dB)', 'FontSize', 12);
end

plot_TFR_across_single(TFR_across, [2:4 8:10 14:16], [5:7 11:13 17:19], 'all', ...
    {"Presentation 1", "Presentation 2", "Presentation 2-1"}, ...
    'TFR across channels and subjects: Comparison between Presentation');
