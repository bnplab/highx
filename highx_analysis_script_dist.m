%% HIGHX Project
% High Excitability State Project
% using post-hoc sorting script from PHASTIMATE code

% Note: Some channels in some subjects have a drift from the preceding TMS
% pulse. We address this using the following measures:
%  - by using a fairly short pre-stimulus window
%  - by applying a linear detrend
%  - by removing drifting channels using a range of mean zscore criterion
% However, this affects the frequency resolution of the spectral
% estimation. An alternative solution would be to use detrending. Note that
% this drift affects the SOUND cleaning, the spectral estimation, but also
% the phase estimation.

% Notes:
% - only uses 64 channels and 800 trials for each dataset

RAW_DATA_FOLDER = 'X:';
TOOLBOXPATH = ['..' filesep 'toolboxes'];

addpath(fullfile(TOOLBOXPATH, 'neurone_tools_for_matlab_1.1.3.11'))
addpath(fullfile(TOOLBOXPATH, 'phastimate'))
addpath(fullfile(TOOLBOXPATH, 'SOUND_functions'))
addpath(fullfile(TOOLBOXPATH, 'circstat-matlab-302acc1'))
addpath(fullfile(TOOLBOXPATH, 'bbci-2e6fe94'))

startup_bbci_toolbox

%% Constants

EEG_EPOCH_LIMITS = [-0.723 -0.005]; %s
EMG_EPOCH_LIMITS = [-0.100  0.100]; %s
EMG_BASELINE_WIN = [-0.100 -0.005]; %s

SPATIAL_FILTER_CHANNELS = {'C3', 'CP1', 'CP5', 'FC1', 'FC5'};
SPATIAL_FILTER_WEIGHTS = [1 -0.25 -0.25 -0.25 -0.25];

%% Initialize Output variables

resultsTable = table();
dataRecordArray = struct([]);

%%

T = readtable('HIGHX_Sessions.xlsx', 'Basic', 1);
assert(~isempty(T), 'table data empty');

T(arrayfun(@(x) strcmp(x, ''), T.sessionFile), :) = []; % remove empty rows

for rowIndex = 1:height(T)

    row = T(rowIndex,:);
    
    fprintf('Processing subject %s ...\n', row.subject{:});

    if (false) % load epoched data
        
        fprintf('Loading epoched data ...')
        load(fullfile(RAW_DATA_FOLDER, '2021-01 HIGHX Epoched Data', '719ms_new', row.subject{:}))
        
    else % generate epoched data
        
        dataDir = fullfile(RAW_DATA_FOLDER, row.project{:}, row.folder{:}, row.subject{:}, row.sessionFile{:}, num2str(row.sessionIndex));
        assert(isfolder(dataDir), 'Data directory not found')

        record = module_read_neurone(convertStringsToChars(fullfile(RAW_DATA_FOLDER, row.project{:}, row.folder{:}, row.subject{:}, row.sessionFile{:}, ['NeurOne-' row.sessionFile{:} '.ses'])), row.sessionIndex);

        fprintf('Pre-processing data ...')

        % Extract timeseries data
        fs = record.properties.samplingRate;
        assert(fs == 5000, 'code is designed to work with EEG data at 5 kHz sampling rate')

        clab = fieldnames(record.signal); %all channel labels
        clab_eeg = clab(structfun(@(chan) strcmp(chan.SignalType, 'EEG'), record.signal));
        clab_emg = clab(structfun(@(chan) strcmp(chan.SignalType, 'EMG'), record.signal));
        data_eeg = zeros(numel(record.signal.(clab_eeg{1}).data), numel(clab_eeg));
        data_emg = zeros(numel(record.signal.(clab_emg{1}).data), numel(clab_emg));
        for i = 1:numel(clab_eeg), data_eeg(:,i) = record.signal.(clab_eeg{i}).data; end
        for i = 1:numel(clab_emg), data_emg(:,i) = record.signal.(clab_emg{i}).data; end
        clear('clab')

        fprintf('.')

        % only consider first 64 channels
        data_eeg(:,65:end) = [];
        clab_eeg(65:end) = [];

        % Create epochs
        ix = record.markers.index(strcmp(record.markers.type, 'Out'));
        ix = ix(ix + min([EEG_EPOCH_LIMITS(1) EMG_EPOCH_LIMITS(1)]) * fs > 0); % keep only indices with sufficient pre-stimulus data
        ix = ix(length(data_eeg) - ix > max([EEG_EPOCH_LIMITS(2) EMG_EPOCH_LIMITS(2)]) * fs); % and with sufficient post-stimulus data
        assert(numel(ix) > 0)
        
        if numel(ix) > 800, ix = ix(1:800); warning('keeping fist 800 trials only'), end;
        
        epochs_eeg = [];
        epochs_eeg.trial = zeros(diff(EEG_EPOCH_LIMITS*fs)+1, length(ix), size(data_eeg, 2));
        epochs_eeg.dimord = 'time_rpt_chan';
        epochs_eeg.time = (EEG_EPOCH_LIMITS(1)*fs:EEG_EPOCH_LIMITS(2)*fs)/fs;
        epochs_eeg.clab = clab_eeg;
        epochs_emg = [];
        epochs_emg.trial = zeros(diff(EMG_EPOCH_LIMITS*fs)+1, length(ix), size(data_emg, 2));
        epochs_emg.dimord = 'time_rpt_chan';
        epochs_emg.time = (EMG_EPOCH_LIMITS(1)*fs:EMG_EPOCH_LIMITS(2)*fs)/fs;
        epochs_emg.clab = clab_emg;
        for e = 1:length(ix) % [time, trial, channel]
            epochs_eeg.trial(:,e,:) = data_eeg(ix(e)+(EEG_EPOCH_LIMITS(1)*fs:EEG_EPOCH_LIMITS(2)*fs), :);
            epochs_emg.trial(:,e,:) = data_emg(ix(e)+(EMG_EPOCH_LIMITS(1)*fs:EMG_EPOCH_LIMITS(2)*fs), :);
        end
        clear('record', 'data_eeg', 'data_emg')

        save(fullfile(RAW_DATA_FOLDER, '2021-01 HIGHX Epoched Data', row.subject{:}), 'epochs_eeg', 'epochs_emg', 'fs', '-v7.3')
    end
    
    % shuffle emg data trials
    %epochs_emg.trial = epochs_emg.trial(:,randperm(size(epochs_emg.trial, 2)), :);

    epochs_eeg.trial = detrend(epochs_eeg.trial, 1); % linear fit de-trend
    epochs_emg.trial = epochs_emg.trial - mean(epochs_emg.trial(epochs_emg.time > EMG_BASELINE_WIN(1) & epochs_emg.time < EMG_BASELINE_WIN(2),:,:), 1); % baseline correction
    
    % Determine MEP amplitudes
    %MEP_LIMITS = cell2num(row.mepLimits)
    MEP_LIMITS = [0.020 0.040];
    epochs_emg.trialinfo = squeeze(range(epochs_emg.trial(epochs_emg.time >= MEP_LIMITS(1) & epochs_emg.time <= MEP_LIMITS(2),:,:)));

    %plot(log(epochs_emg.trialinfo(:,1)), log(epochs_emg.trialinfo(:,2)), '.')
    %loglog(epochs_emg.trialinfo(:,1), epochs_emg.trialinfo(:,2), 'o') 
    
    fprintf('.')
    
    % De-mean
    epochs_eeg.trial = detrend(epochs_eeg.trial, 0);
    
    fprintf(' Done.')
    
    % Reject channels with an average range indicating slow drifts that
    % could not be removed with linear detrending
    
    fprintf('\nBad channel and bad trial detection...\n')
    
    %TODO: use range of mean trial in addition to median of trial ranges?
    r = squeeze(median(range(epochs_eeg.trial, 1), 2));
    badchan_logical = r > 150;
    badchan_clab = epochs_eeg.clab(badchan_logical);
    if any(badchan_logical)
        figure;
        h = bar(reordercats(categorical(epochs_eeg.clab), epochs_eeg.clab), r);
        h.FaceColor = 'flat'; h.CData(badchan_logical,1) = 1;
        %title(h, row.subject{:}, 'Interpreter', 'none')
        %ylabel(h, 'median range (µV)')
    end
    if any(badchan_logical)
        warning(['Removing channels due to range criterion:' sprintf(' %s', epochs_eeg.clab{badchan_logical})])        
        epochs_eeg.trial(:,:,badchan_logical) = [];
        epochs_eeg.clab(badchan_logical) = [];
    end
    if ~all(ismember(SPATIAL_FILTER_CHANNELS, epochs_eeg.clab))
        warning('Channels required for spatial filter were removed, skipping this dataset')
        continue
    end
    
    % Bad trial removal
    r = squeeze(max(range(epochs_eeg.trial, 1), [], 3));
    badtrial_logical = r > 500;
    
    fprintf('Removing %i/%i (%.2f%%) bad trials due to range criterion', sum(badtrial_logical == 1), numel(badtrial_logical), 100*sum(badtrial_logical == 1)/numel(badtrial_logical))
    epochs_eeg.trial(:,badchan_logical,:) = [];
    epochs_emg.trial(:,badchan_logical,:) = [];
    epochs_emg.trialinfo(badchan_logical,:) = [];
    
    fprintf('\nDone.')
    
    % convert MEP size to first component of PCA
    % note that this data is then mean-centered (TODO: think about making
    % it median centered)
    [~, score] = pca(log(epochs_emg.trialinfo)); % channels are columns
    excitability_index = score(:,1); % use the first PCA score
    
    fprintf('\nDown-sampling EEG data...')
    % Down-sample eeg data
    % TODO: Make sure that the signal processing filtfilt function is used and not a fieldtrip replacement
    D = designfilt('lowpassfir', 'FilterOrder', 50, 'CutoffFrequency', 250, 'SampleRate', fs);
    epochs_eeg.trial = filtfilt(D, epochs_eeg.trial);
    epochs_eeg.trial = epochs_eeg.trial(1:5:end,:,:);
    epochs_eeg.time = epochs_eeg.time(1:5:end);
    
    fprintf(' Done.\n')
    
    % Prepare leadfield matrix    
    [x,y,z,clab]= mntutil_posExt55;
    [lia, ~] = ismember(epochs_eeg.clab, clab);
    if ~all(lia), warning('removing channels not found in leadfield'), end
    epochs_eeg.trial(:, :, lia == 0) = []; % remove channels that are not in the leadfield
    epochs_eeg.clab(lia == 0) = [];
	[~, locb] = ismember([epochs_eeg.clab; 'FCz'], clab); % location of the reference channel is needed
    LFM = ComputeSphericalLFM_chanlocs_simple(x(locb)', y(locb)', z(locb)', numel(locb)); % lead-field matrix (in same reference as data)
    
    fprintf('Cleaning data with SOUND algorithm ...')
    [corrected_data, sigmas, dn, correctionM] = SOUND_fast_correction(permute(epochs_eeg.trial, [3 1 2]), LFM, 10, 0.01, []); % data: channels x time points x trials (note: use single-channel reference, _not_ average reference)
    epochs_eeg.trial_sound = permute(corrected_data, [2 3 1]);
    fprintf(' Done.\n')
    clear('corrected_data');    
    
    % Create spatial filter (C3 Hjorth)
    assert(all(ismember(SPATIAL_FILTER_CHANNELS, epochs_eeg.clab)), 'not all channels from spatial filter found')
    spatial_filter = zeros(numel(epochs_eeg.clab), 1);
    spatial_filter(ismember(epochs_eeg.clab, SPATIAL_FILTER_CHANNELS)) = SPATIAL_FILTER_WEIGHTS;

    [a, aSOUND, wBF, wSOUND] = getFilterBFPattern_SOUND(permute(epochs_eeg.trial, [3 1 2]), LFM, spatial_filter, correctionM, 0.01, 0.01);
    
    %% Process different spatial filters
    
    % set of spatial filters to apply (to uncleaned data)
    spatial_filter_conditions = {};
    spatial_filter_conditions(1).name = 'RAW_C3LAP';
    spatial_filter_conditions(1).w = spatial_filter;
    spatial_filter_conditions(2).name = 'SOUND_C3LAP';
    spatial_filter_conditions(2).w = correctionM * spatial_filter;
    spatial_filter_conditions(3).name = 'RAW_C3BF';
    spatial_filter_conditions(3).w = wBF;
    spatial_filter_conditions(4).name = 'SOUND_C3BF';
    spatial_filter_conditions(4).w = correctionM * wBF;

    legendNames = {'RAW⋅SL', 'SOUND⋅SL', 'RAW⋅BF', 'SOUND⋅BF'};

    % scale each filter
    for i = 1:numel(spatial_filter_conditions)
        spatial_filter_conditions(i).w = spatial_filter_conditions(i).w ./ norm(spatial_filter_conditions(i).w, 'fro');
    end
    
    epochs_spf = struct();
    freq_spf = struct();
    mdlfit_spf = struct();

    for condition_ix = 1:numel(spatial_filter_conditions)
    
        % Apply spatial filter
        epochs_spf(condition_ix).trial = reshape(epochs_eeg.trial, [], size(epochs_eeg.trial, 3)) * spatial_filter_conditions(condition_ix).w;
        epochs_spf(condition_ix).trial = reshape(epochs_spf(condition_ix).trial, size(epochs_eeg.trial, 1), size(epochs_eeg.trial, 2), []);
        epochs_spf(condition_ix).dimord = 'time_rpt_chan';
        epochs_spf(condition_ix).time = epochs_eeg.time;
        epochs_spf(condition_ix).label = spatial_filter_conditions(condition_ix).name;

        epochs_spf(condition_ix).trial = detrend(epochs_spf(condition_ix).trial, 0); % de-mean (also applied to the shorter pre-stimulus window)

        % NOTE: don't remove bad trials (as different trials could be
        % removed for different filters)

        %% Phase Estimation
        % parameters differ from https://doi.org/10.1016/j.neuroimage.2020.116761 to avoid inhomogenous distribution of phase estimates
        D = designfilt('bandpassfir', 'FilterOrder', 192, 'CutoffFrequency1', 9, 'CutoffFrequency2', 13, 'SampleRate', 1000, 'DesignMethod', 'window');
        epochs_spf(condition_ix).estphase = phastimate(detrend(double(epochs_spf(condition_ix).trial(1:end,:)),0), D, 65, 25, 128, 4);

        %wrap to -1.5*pi .. 0.5*pi (for visualization later)
        epochs_spf(condition_ix).estphase(epochs_spf(condition_ix).estphase > 0.5*pi) = epochs_spf(condition_ix).estphase(epochs_spf(condition_ix).estphase > 0.5*pi)-(2*pi);
 
        % Create trialmask variable, allowing exclusion of trials
        trialmask = true(1, size(excitability_index, 1));

        % Apply percentile filter for logistic regression
        exc_low = excitability_index(trialmask) - RankOrderFilter(excitability_index(trialmask), 25, 25) < 0;
        exc_high = excitability_index(trialmask) - RankOrderFilter(excitability_index(trialmask), 25, 75) > 0;    
        y_binary = nan(size(excitability_index));
        trialmask_binary = trialmask;
        trialmask_binary(trialmask) = exc_high | exc_low;
        y_binary(trialmask) = exc_high - exc_low;

        mdl_linear.name = 'linear';
        mdl_linear.estphase = epochs_spf(condition_ix).estphase(trialmask);
        mdl_linear.excitability = excitability_index(trialmask)';
        mdl_linear.fit = fitlm([cos(epochs_spf(condition_ix).estphase(trialmask))' sin(epochs_spf(condition_ix).estphase(trialmask))'], excitability_index(trialmask));
        mdl_linear.shuffled = {};
        for rep = 1:100
            y_shuffled = excitability_index(trialmask);
            y_shuffled = y_shuffled(randperm(size(y_shuffled, 1)));
            mdl_linear.shuffled{rep} = fitlm([cos(epochs_spf(condition_ix).estphase(trialmask))' sin(epochs_spf(condition_ix).estphase(trialmask))'], y_shuffled);
        end
        
        mdl_twoclass.name = 'twoclass';
        mdl_twoclass.estphase = epochs_spf(condition_ix).estphase(trialmask_binary);
        mdl_twoclass.excitability = y_binary(trialmask_binary)';
        mdl_twoclass.fit = fitlm([cos(epochs_spf(condition_ix).estphase(trialmask_binary))' sin(epochs_spf(condition_ix).estphase(trialmask_binary))'], y_binary(trialmask_binary));
        mdl_twoclass.shuffled = {};
        for rep = 1:100
            y_shuffled = y_binary(trialmask_binary);
            y_shuffled = y_shuffled(randperm(size(y_shuffled, 1)));
            mdl_twoclass.shuffled{rep} = fitlm([cos(epochs_spf(condition_ix).estphase(trialmask_binary))' sin(epochs_spf(condition_ix).estphase(trialmask_binary))'], y_shuffled);
        end
        
        mdlfit_spf(condition_ix).label = epochs_spf(condition_ix).label;
        mdlfit_spf(condition_ix).weights = {spatial_filter_conditions(condition_ix).w};
        mdlfit_spf(condition_ix).meantrial = {mean(epochs_spf(1).trial, 2)};
        mdlfit_spf(condition_ix).mediantrial = {median(epochs_spf(1).trial, 2)};
        mdlfit_spf(condition_ix).estphase = {epochs_spf(condition_ix).estphase}; % phase estimates for all trials
        mdlfit_spf(condition_ix).snr = {freq_spf(condition_ix).snr};
        
        % Test for homogeneity
        mdlfit_spf(condition_ix).rtest = circ_rtest(epochs_spf(condition_ix).estphase);
            
        for mdl = {mdl_linear, mdl_twoclass}
        
            coeffs = mdl{:}.fit.Coefficients.Estimate;

            % determine the phase at which the excitability is hightest:
            % to understand the below, see http://scipp.ucsc.edu/~haber/ph5B/addsine.pdf
            phi = wrapToPi(sign(sin(coeffs(3)/sqrt(coeffs(2)^2 + coeffs(3)^2))) * pi/2 - atan(coeffs(2)/coeffs(3)));

            %wrap to -1.5*pi .. 0.5*pi
            if phi > 0.5*pi, phi = phi-(2*pi); end

    %         bestfit = @(x) coeffs(1) + coeffs(2)*cos(x) + coeffs(3)*sin(x);
    %         figure, plot((-1.5*pi):0.01:(0.5*pi), bestfit((-1.5*pi):0.01:(0.5*pi)))
    %         xline(phi)     

            fprintf('%s fit has p-value = %.3f and peak excitability at phi = %3.0f° (%s)\n', mdl{:}.name, mdl{:}.fit.coefTest, rad2deg(phi), epochs_spf(condition_ix).label)
            
            mdlfit_spf(condition_ix).([mdl{:}.name '_coeff1']) = coeffs(1);
            mdlfit_spf(condition_ix).([mdl{:}.name '_coeff2']) = coeffs(2);
            mdlfit_spf(condition_ix).([mdl{:}.name '_coeff3']) = coeffs(3);
            mdlfit_spf(condition_ix).([mdl{:}.name '_pValue']) = mdl{:}.fit.coefTest;
            mdlfit_spf(condition_ix).([mdl{:}.name '_phi_deg']) = rad2deg(phi);
            mdlfit_spf(condition_ix).([mdl{:}.name '_Rsquared_ord']) = mdl{:}.fit.Rsquared.Ordinary;
            mdlfit_spf(condition_ix).([mdl{:}.name '_Rsquared_adj']) = mdl{:}.fit.Rsquared.Adjusted;
            
            coeffs_shuffled1 = cellfun(@(x) x.Coefficients.Estimate(1), mdl{:}.shuffled, 'UniformOutput', true);
            coeffs_shuffled2 = cellfun(@(x) x.Coefficients.Estimate(2), mdl{:}.shuffled, 'UniformOutput', true);
            coeffs_shuffled3 = cellfun(@(x) x.Coefficients.Estimate(3), mdl{:}.shuffled, 'UniformOutput', true);
            phi_shuffled = wrapToPi(sign(sin(coeffs_shuffled3./sqrt(coeffs_shuffled2.^2 + coeffs_shuffled3.^2))) .* pi/2 - atan(coeffs_shuffled2./coeffs_shuffled3));
            
            mdlfit_spf(condition_ix).([mdl{:}.name '_shuffled_coeff1']) = coeffs_shuffled1;
            mdlfit_spf(condition_ix).([mdl{:}.name '_shuffled_coeff2']) = coeffs_shuffled2;
            mdlfit_spf(condition_ix).([mdl{:}.name '_shuffled_coeff3']) = coeffs_shuffled3;
            mdlfit_spf(condition_ix).([mdl{:}.name '_shuffled_pValue']) = cellfun(@(x) x.coefTest, mdl{:}.shuffled, 'UniformOutput', true);
            mdlfit_spf(condition_ix).([mdl{:}.name '_shuffled_phi_deg']) = rad2deg(phi_shuffled);
            mdlfit_spf(condition_ix).([mdl{:}.name '_shuffled_Rsquared_ord']) = cellfun(@(x) x.Rsquared.Ordinary, mdl{:}.shuffled, 'UniformOutput', true);
            mdlfit_spf(condition_ix).([mdl{:}.name '_shuffled_Rsquared_adj']) = cellfun(@(x) x.Rsquared.Adjusted, mdl{:}.shuffled, 'UniformOutput', true);
            
        end % model

    end % spatial filter condition

    %% append data to table
     
    % append model data for each spatial filter
    for condition_ix = 1:numel(mdlfit_spf)
        row = [row renamevars(struct2table(mdlfit_spf(condition_ix)), 1:width(struct2table(mdlfit_spf(condition_ix))), strcat(mdlfit_spf(condition_ix).label, '_', struct2table(mdlfit_spf(condition_ix)).Properties.VariableNames))];
    end
    
    row.numTrials = size(epochs_eeg.trial, 2);
    row.badchan_num = sum(badchan_logical == 1);
    row.badchan_clab = {badchan_clab};
    %row.badtrial_num = sum(badtrial_logical == 1); % currently there is no bad trial rejection
    
    row.clab = {epochs_eeg.clab};
    row.spatial_filters = spatial_filter_conditions;

    row.freqboi_cf = freqboi_cf;

    %% append variables to struct
    
    resultsTable = [resultsTable; row];
    
    
end

save resultsTable
save(['resultsTable_' datestr(now, 'yyyy-mm-dd') '.mat'], 'resultsTable');

copyfile(getfield(matlab.desktop.editor.getActive, 'Filename'), [fileparts(getfield(matlab.desktop.editor.getActive, 'Filename')) filesep 'highx_mainscript_datarun_' datestr(now, 'yyyy-mm-dd') '.m']);

writetable(resultsTable, ['resultsTable_' datestr(now, 'yyyy-mm-dd') '.xlsx']);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ANALYSIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

set(0,'DefaultFigureWindowStyle','normal') 

ALPHABET = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ';

mdlNames = {'RAW_C3LAP', 'SOUND_C3LAP', 'RAW_C3BF', 'SOUND_C3BF'};
legendNames = {'RAW⋅SL', 'SOUND⋅SL', 'RAW⋅BF', 'SOUND⋅BF'};


%% calculations
for i = 1:numel(mdlNames)
    spect_snr = cat(2, resultsTable.([mdlNames{i} '_snr']){:});
    freqboi_cf = resultsTable(1,:).freqboi_cf;
    
    alpha_bin_indices = find(freqboi_cf >= 8.5 & freqboi_cf <= 14);
    % detect peak frequencies
    diff_spect_snr = diff(spect_snr); % run diff with the whole data to avoid edge effects in region of interest
    has_alpha_peak = true(1, size(spect_snr, 2));
    has_alpha_peak = has_alpha_peak & (diff_spect_snr(alpha_bin_indices(1),:) > 0 & diff_spect_snr(alpha_bin_indices(end),:) < 0); % rising at beginning, falling at end
    has_alpha_peak = has_alpha_peak & (sum(diff(sign(diff_spect_snr(alpha_bin_indices,:))) == -2) == 1); % has single peak
    
    mu_snr_db = nan(size(spect_snr, 2), 1);
    mu_iaf = nan(size(spect_snr, 2), 1);
    [snr, ix] = max(spect_snr(alpha_bin_indices, has_alpha_peak));
    mu_snr_db(has_alpha_peak) = pow2db(snr);
    mu_iaf(has_alpha_peak) = freqboi_cf(alpha_bin_indices(ix));
    
    snrTable = table(mu_iaf, mu_snr_db);
    snrTable.Properties.VariableNames = strcat([mdlNames{i}, '_'], snrTable.Properties.VariableNames);
    
    resultsTable = [resultsTable snrTable];
end


%% Methodological figure about artifacts in hilbert transform

f = figure('Color', 'white');
f.Units = 'centimeters';
f.Position = [0 0 24 15];

t = -0.5:0.001:0.5;
x = sin( t .* 3 .* (2*pi) )';
x = x + [0 1 -1] * 0.5; % add offset
phi = angle(hilbert(x));

ax = subplot(2, 3, [1:3]);
plot(ax, t, x, 'LineWidth', 2); hold on
ax.ColorOrderIndex = 1;
plot(ax, t, angle(hilbert(x)), 'LineWidth', 2, 'LineStyle', ':')
yline(ax, 0, 'LineWidth', 1, 'LineStyle', '--')
axis(ax,'tight')
title(ax, 'A'), ax.TitleHorizontalAlignment = 'left';

h1 = subplot(2, 3, 4, polaraxes);
polarhistogram(h1, phi(:, 1), 12, 'FaceColor', ax.ColorOrder(1,:))
title(h1, 'B'), h1.TitleHorizontalAlignment = 'left';

h2 = subplot(2, 3, 5, polaraxes);
polarhistogram(h2, phi(:, 2), 12, 'FaceColor', ax.ColorOrder(2,:))
title(h2, 'C'), h2.TitleHorizontalAlignment = 'left';

h3 = subplot(2, 3, 6, polaraxes);
polarhistogram(h3, phi(:, 3), 12, 'FaceColor', ax.ColorOrder(3,:))
title(h3, 'D'), h3.TitleHorizontalAlignment = 'left';

for h = [h1 h2 h3]
    h.ThetaAxisUnits = 'radians';
    h.ThetaLim = [-pi pi];
    h.ThetaZeroLocation = 'top';
    h.ThetaDir = 'clockwise';
end

%% Another one with a full-wave rectified sine

f = figure('Color', 'white');
f.Units = 'centimeters';
f.Position = [0 0 24 15];

t = -0.5:0.001:0.5;
x = abs(sin( t .* 3/2 .* (2*pi) ))';
x = x - mean(x)
x = x + [0 1 -1] * 0.25; % add offset
phi = angle(hilbert(x));

ax = subplot(2, 3, [1:3]);
plot(ax, t, x, 'LineWidth', 2); hold on
ax.ColorOrderIndex = 1;
plot(ax, t, angle(hilbert(x)), 'LineWidth', 2, 'LineStyle', ':')
yline(ax, 0, 'LineWidth', 1, 'LineStyle', '--')
axis(ax,'tight')
title(ax, 'A'), ax.TitleHorizontalAlignment = 'left';

h1 = subplot(2, 3, 4, polaraxes);
polarhistogram(h1, phi(:, 1), 24, 'FaceColor', ax.ColorOrder(1,:))
title(h1, 'B'), h1.TitleHorizontalAlignment = 'left';

h2 = subplot(2, 3, 5, polaraxes);
polarhistogram(h2, phi(:, 2), 24, 'FaceColor', ax.ColorOrder(2,:))
title(h2, 'C'), h2.TitleHorizontalAlignment = 'left';

h3 = subplot(2, 3, 6, polaraxes);
polarhistogram(h3, phi(:, 3), 24, 'FaceColor', ax.ColorOrder(3,:))
title(h3, 'D'), h3.TitleHorizontalAlignment = 'left';

for h = [h1 h2 h3]
    h.ThetaAxisUnits = 'radians';
    h.ThetaLim = [-pi pi];
    h.ThetaZeroLocation = 'top';
    h.ThetaDir = 'clockwise';
end


%% Check for homogeneity and plot overall result

EDGES = -pi:pi/12:pi;

%linearOrTwoclass = 'twoclass';
linearOrTwoclass = 'linear';

f = figure('Color', 'white');
f.Units = 'centimeters';
f.Position = [0 0 23 15];

ax1 = subplot(2,3,[1 2 4 5],polaraxes); hold(ax1, 'on');
ax2 = subplot(2,3,3,polaraxes); hold(ax2, 'on');
ax3 = subplot(2,3,6,polaraxes); hold(ax3, 'on');
for mdlName = mdlNames(1)
    includeMask = true(height(resultsTable), 1);
    %includeMask = includeMask & includeMask_rtest & includeMask_pvalue;
    
    estphase_data = resultsTable.([mdlName{:} '_estphase'])(includeMask);
    estphase_data = cat(2,  estphase_data{:})';
    shuffled_phi_data = deg2rad(resultsTable.([mdlName{:} '_' linearOrTwoclass '_shuffled_phi_deg'])(includeMask, :));
    shuffled_phi_data = shuffled_phi_data(:);
    phi_data = deg2rad(resultsTable.([mdlName{:} '_' linearOrTwoclass '_phi_deg'])(includeMask));

    fprintf('\n%s:\n', mdlName{:});
    fprintf('estphase rtest p-value = %.3f\n', circ_rtest(estphase_data));
    fprintf('shuffled rtest p-value = %.3f\n', circ_rtest(shuffled_phi_data));
    fprintf('phi_data rtest p-value = %.1e\n', circ_rtest(phi_data));
    fprintf('phi_data mean %.2f (median %.2f) differs from pi at alpha 1e-3 mtest rejected = %i\n', rad2deg(circ_mean(phi_data)), rad2deg(circ_median(phi_data)), circ_mtest(phi_data, pi, 1e-3));
    
    polarhistogram(ax1, phi_data, EDGES, 'Normalization', 'count')%, 'DisplayStyle', 'stairs', 'LineWidth', 1.5)
    polarhistogram(ax2, estphase_data, EDGES, 'Normalization', 'pdf')%, 'DisplayStyle', 'stairs', 'LineWidth', 1.5)
    polarhistogram(ax3, shuffled_phi_data, EDGES, 'Normalization', 'pdf')%, 'DisplayStyle', 'stairs', 'LineWidth', 1.5)
        
    rlim_max = ax1.RLim(2) + 1;
    text(ax1, pi/2, rlim_max/2, sprintf('n=%i', numel(phi_data)));    
    
    [alpha_bar ul ll] = circ_mean(phi_data);
    [alpha_bar] = circ_median(phi_data);
    polarplot(ax1, [1 1] .* alpha_bar, [0 1] .* rlim_max, 'r-', 'LineWidth', 2) % median
    polarplot(ax1, linspace(ll, ul, 20), repmat(1, 1, 20) .* rlim_max, 'r-', 'LineWidth', 2) % variance
    
    polarplot(ax1, -pi + 2*pi*11.5*[0.000 0.010 0.020], ax1.RAxis.Limits(2) * 1.00, 'ko', 'MarkerSize', 7.5, 'LineWidth', 1) % 3 100 Hz stimuli assuming 11.5 Hz peak frequency
    rad2deg(ll)
    rad2deg(ul)
end
title(ax1, 'A'), ax1.TitleHorizontalAlignment = 'left';
title(ax2, 'B'), ax2.TitleHorizontalAlignment = 'left';
title(ax3, 'C'), ax3.TitleHorizontalAlignment = 'left';

text(ax1, pi/4, ax1.RLim(2)*.7, sprintf('n = %i', sum(includeMask)), 'HorizontalAlignment', 'center', 'FontSize', 12);

%l = legend(ax1, legendNames, 'Interpreter', 'none', 'Location', 'east');

% equalize scales, rotate, change labels
for h = [ax1 ax2 ax3]
    h.ThetaAxisUnits = 'radians';
    h.ThetaLim = [-pi pi];
    h.ThetaZeroLocation = 'top';
    h.ThetaDir = 'clockwise';
end

ax1.ThetaTick = (-5/6*pi):pi/6:pi;
ax1.ThetaTickLabel = {'-5\pi/6', '-2\pi/3', 'rising', '-\pi/3', '-\pi/6', 'positive peak', '\pi/6', '\pi/3', 'falling', '2\pi/3', '5\pi/6', 'negative peak'};

ax2.ThetaTickLabel = {};
ax3.ThetaTickLabel = {};


%% Visualize Group Results

%linearOrTwoclass = 'twoclass';
linearOrTwoclass = 'linear';

% inclusion criteria:
fprintf('\nExcluding subjects:')
% (1) with an inhomogenous phase distribution
includeMask_rtest = true(height(resultsTable), 1);
for mdlName = mdlNames
    includeMask_rtest = includeMask_rtest & resultsTable.([mdlName{:} '_rtest']) > 0.05;
end
% (2) without non-significant phase effect in all filters
includeMask_pvalue = false(height(resultsTable), 1); % start excluding
for mdlName = mdlNames
    includeMask_pvalue = includeMask_pvalue | resultsTable.([mdlName{:} '_' linearOrTwoclass '_pValue']) < 0.05;
end
fprintf(' %i of %i subjects excluded due to inhomogenous phase distribution (%i) or insignificant phase effect (%i)\n', ...
    sum((includeMask_rtest & includeMask_pvalue) == 0), ...
    numel(includeMask_rtest & includeMask_pvalue), ...
    sum(includeMask_rtest == 0), ...
    sum(includeMask_pvalue == 0));

% Peak excitability phase plot

fig = figure('Color', 'white');
fig.Units = 'centimeters';
fig.Position = [0 0 40 18];

EDGES = -pi:pi/16:pi;

h1 = {}; h2 = {}; t = {};

for i = 1:numel(mdlNames)
    h1{i} = subplot(2, numel(mdlNames),                 i, polaraxes, 'Parent', fig); hold on
    h2{i} = subplot(2, numel(mdlNames), numel(mdlNames)+i, polaraxes, 'Parent', fig); hold on
    
    % option (1): exclude by subject
    h1{i}.UserData = deg2rad(resultsTable.([mdlNames{i} '_' linearOrTwoclass '_phi_deg'])( includeMask_rtest & includeMask_pvalue ));
    % option (2): exclude individual filters
    h2{i}.UserData = deg2rad(resultsTable.([mdlNames{i} '_' linearOrTwoclass '_phi_deg'])( (resultsTable.([mdlNames{i} '_rtest']) > 0.05) & (resultsTable.([mdlNames{i} '_' linearOrTwoclass '_pValue']) < 0.05) ));

    fprintf('\n%s by-subject   mtest that mean phase %.2f distribution contains pi at alpha 0.05 rejected = %i', mdlNames{i}, rad2deg(circ_mean(h1{i}.UserData)), circ_mtest(h1{i}.UserData, pi, 0.05));
    fprintf('\n%s by-transform mtest that mean phase %.2f distribution contains pi at alpha 0.05 rejected = %i', mdlNames{i}, rad2deg(circ_mean(h2{i}.UserData)), circ_mtest(h2{i}.UserData, pi, 0.05));

    polarhistogram(h1{i}, h1{i}.UserData, EDGES, 'FaceColor', h1{i}.ColorOrder(i,:))
    polarhistogram(h2{i}, h2{i}.UserData, EDGES, 'FaceColor', h2{i}.ColorOrder(i,:))
    
end

fig.Units = 'normalized';
for i=1:numel(mdlNames)
    annotation('textbox', [h1{i}.Position(1), 0.475, h1{i}.Position(3)-0.02, 0.05], 'string', legendNames{i}, 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'LineStyle', 'none')
end
annotation('textbox', [0, h1{1}.Position(2), h1{1}.OuterPosition(1), h1{1}.Position(4)], 'string', 'A', 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'LineStyle', 'none')
annotation('textbox', [0, h2{1}.Position(2), h2{1}.OuterPosition(1), h2{1}.Position(4)], 'string', 'B', 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'LineStyle', 'none')

% average mean angle
cellfun(@(x) rad2deg(circ_mean(x.UserData)), [h1 h2])
rad2deg(circ_mean(cellfun(@(x) circ_mean(x.UserData), [h1 h2])))

% average median angle
cellfun(@(x) rad2deg(circ_median(x.UserData)), [h1 h2])
rad2deg(circ_mean(cellfun(@(x) circ_median(x.UserData), [h1 h2])))

for h = [h1 h2]
    
    h{:}.ThetaAxisUnits = 'radians';
    h{:}.ThetaLim = [-pi pi];
    h{:}.ThetaZeroLocation = 'top';
    h{:}.ThetaDir = 'clockwise';

    h{:}.ThetaTick = pi .* [(-6+1):2*6] ./ 6;
    h{:}.ThetaTickLabel = repmat({''}, 1, numel(h{:}.ThetaTick));
    
    %labels only on the outside
    if h{:}.OuterPosition(2) > 0.5, h{:}.ThetaTickLabel(h{:}.ThetaTick == -0) = {'pos peak'}; end % top row
    if h{:}.OuterPosition(2) < 0.5, h{:}.ThetaTickLabel(h{:}.ThetaTick == pi) = {'neg peak'}; end % bottom row
    if h{:}.OuterPosition(1) < 1/numel(mdlNames), h{:}.ThetaTickLabel(h{:}.ThetaTick == -pi/2) = {'rising'}; end % left column
    if (h{:}.OuterPosition(1) + h{:}.OuterPosition(3))  > 1-1/numel(mdlNames), h{:}.ThetaTickLabel(h{:}.ThetaTick == pi/2) = {'falling'}; end % left column
    
    rlim_max = h{:}.RLim(2) + 1;
    text(h{:}, pi/2, rlim_max/2, sprintf('n=%i', numel(h{:}.UserData)));
    
    [alpha_bar] = circ_median(h{:}.UserData);
    [~, ul, ll] = circ_mean(h{:}.UserData);
    polarplot(h{:}, [1 1] .* alpha_bar, [0 1] .* rlim_max, 'r-', 'LineWidth', 2) % median
    polarplot(h{:}, linspace(ll, ul, 20), repmat(1, 1, 20) .* rlim_max, 'r-', 'LineWidth', 2) % variance
    
    polarplot(h{:}, -pi + 2*pi*11.5*[0.000 0.010 0.020], h{:}.RAxis.Limits(2) * 1.00, 'ko', 'MarkerSize', 7.5, 'LineWidth', 1) % 3 100 Hz stimuli assuming 11.5 Hz peak frequency
    %polarplot(h{:}, -pi + 2*pi*13*[0.000 0.010 0.020], h{:}.RAxis.Limits(2)* 1.00, 'k*', 'MarkerSize', 7.5, 'LineWidth', 1) % 3 100 Hz stimuli assuming 13 Hz peak frequency
end  

% Circular statistics

pair_ix = nchoosek(1:numel(mdlNames), 2);

pvalue_data = nan(numel(mdlNames));
for i = 1:size(pair_ix, 1)
    [pvalue_data(pair_ix(i, 2), pair_ix(i, 1)) ~] = circ_wwtest(h2{pair_ix(i, 1)}.UserData, h2{pair_ix(i, 2)}.UserData);
end

t = array2table(pvalue_data);
t.Properties.VariableNames = mdlNames;
t.Properties.RowNames = mdlNames;
t.Properties.Description = 'Parametric Watson-Williams multi-sample test for equal means';

% Compare filters by pValue and Rsquared

fig = figure('Color', 'white');
%fig.WindowStyle = 'docked';
fig.Units = 'centimeters';
fig.Position = [0 0 20 10];


includeMask = true(height(resultsTable), 1);
includeMask = includeMask & includeMask_rtest; % only homogeneous distributions

ax = subplot(1,2,1); hold on

for mdlName = strcat(mdlNames, '_', linearOrTwoclass)
    histogram(resultsTable.([mdlName{:} '_pValue'])(includeMask), 'Normalization', 'cdf', 'BinWidth', 0.00001, 'DisplayStyle', 'stairs', 'LineWidth', 1.5)
end
legend(legendNames, 'Interpreter', 'none', 'Location', 'southeast', 'AutoUpdate', 'off')

ax.YLim = [0 0.75];
ax.XLim = [0.0000 0.5];
%ax.XScale = 'log';
ax.XTick = [ 0.1 0.2 0.3 0.4 0.5];
ax.YTickLabel = 100 * ax.YTick;
ax.XLabel.String = 'p-value';
ax.YLabel.String = 'Participant Percentile';

yline(0.5)
ax.ColorOrderIndex = 1;
for mdlName = strcat(mdlNames, '_', linearOrTwoclass)
    plot(ax, median(resultsTable.([mdlName{:} '_pValue'])(includeMask)), 0.5, '.', 'MarkerSize', 25)
end

xline(0.1)
ax.ColorOrderIndex = 1;
fprintf('Percentile reaching a p-Value of 0.1')
for mdlName = strcat(mdlNames, '_', linearOrTwoclass)
    proportion = sum(resultsTable.([mdlName{:} '_pValue'])(includeMask) < 0.1) / sum(includeMask);
    fprintf('%s: %.2f\n', mdlName{:}, proportion)  
    plot(ax, 0.1, proportion, '.', 'MarkerSize', 25)
end

%title(['median p-value in ', linearOrTwoclass, ' model'])
title('A')
ax.TitleHorizontalAlignment = 'left';

ax = subplot(1,2,2); hold on

for mdlName = strcat(mdlNames, '_', linearOrTwoclass)
    histogram(resultsTable.([mdlName{:} '_Rsquared_ord'])(includeMask), 'Normalization', 'cdf', 'BinWidth', 0.00001, 'DisplayStyle', 'stairs', 'LineWidth', 1.5)
end

ax.YLim = [0 0.75];
ax.XLim = [0.000 .005];
%ax.XScale = 'log';
ax.XTick = [0.000:0.001:0.005];
ax.YTickLabel = 100 * ax.YTick;
ax.XLabel.String = 'R^2';

yline(0.5)

ax.ColorOrderIndex = 1;
for mdlName = strcat(mdlNames, '_', linearOrTwoclass)
    plot(ax, median(resultsTable.([mdlName{:} '_Rsquared_ord'])(includeMask)), 0.5, '.', 'MarkerSize', 25)
end

%title(['median Rsquared in ', linearOrTwoclass, ' model'])
title('B')
ax.TitleHorizontalAlignment = 'left';

%% Spatial Filter and Signal Analysis
%TODO: refactor following spatial filter name as column name convention

figure;
ax = {};

for i = 1:numel(mdlNames)

    spect_snr = cat(2, resultsTable.([mdlNames{i} '_snr']){:});
    freqboi_cf = resultsTable(1,:).freqboi_cf;
    
    %spect_snr = spect_snr(:, has_alpha_peak);
    %snr = snr(:, includeMask_rtest & includeMask_pvalue);

    ax{i} = subplot(1, numel(mdlNames), i);
    hold(ax{i}, 'on');
    data = pow2db(spect_snr);
    plot(ax{i}, freqboi_cf, data, 'Color', [1 1 1] .* 0.67)
    plot(ax{i}, resultsTable.([mdlNames{i}, '_mu_iaf']), resultsTable.([mdlNames{i}, '_mu_snr_db']), 'r*')
    plot(ax{i}, freqboi_cf, mean(data, 2), 'LineWidth', 2, 'Color', 'k')
    ax{i}.XScale = 'log';
    ax{i}.XTick = [5 10 20 50];
    xlabel(ax{i}, 'Hz')
    ax{i}.XLim = [5 60];
    ax{i}.YLim = [-1 20];
    title(ax{i}, legendNames{i})
    
end
linkaxes([ax{:}], 'y')
ylabel(ax{1}, 'SNR (dB)')
