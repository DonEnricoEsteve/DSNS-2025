%% load data (averaged across trials per condition per subject)

load("don_freq_across_allcons.mat")
load("don_freq_across_allconsbaseline.mat")

%% Normalize and extract powspctrm values from frequency data 
function avgPower = extract_2D_power(baseline_cell, poststim_cell, channels, freq_ranges)
    
    [nSubj, nCond] = size(poststim_cell);
    avgPower = cell(nSubj, nCond);

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
            

            % --- Average within the frequency bands ---
            mat_av_bands = [];

            for freqrange=freq_ranges
              
                fidx = postSel.freq >= freqrange{1,1}(1) & postSel.freq <= freqrange{1,1}(2);
                mat_av_bands = [mat_av_bands, mean(normPow(:,fidx), 2, 'omitnan')];

            end

            avgPower(s, c) = {mat_av_bands};
        end
    end
end

tic

freq_ranges = {[4 8], [8 13], [13 25], [25 50]};
avg2DPower = extract_2D_power(freqbase_across, freq_across, "all", freq_ranges);


toc

% Global (all channels)
% θ-band (4–8 Hz) 
% α-band (8–13 Hz) 
% β-band (13–25 Hz)
% γ-band (25–50 Hz)

avg2DPower_struct.data = avg2DPower;
avg2DPower_struct.freqs = freq_across{1,1}.freq;
avg2DPower_struct.chan_labels = freq_across{1,1}.label;

save("avg2DPower.mat", "avg2DPower_struct")


%% Across Subject Decoding %%

%% Preparing data %%

% load the avg2DPower structure
load("avg2DPower.mat")

% extract class labels and data
data = avg2DPower_struct.data;

u_labels = categorical({'odd', 'FS1', 'FM1', 'FL1', 'FS2', 'FM2', 'FL2', ...
    'PS1', 'PM1', 'PL1', 'PS2', 'PM2', 'PL2', ...
    'NS1', 'NM1', 'NL1', 'NS2', 'NM2', 'NL2'});

% reshape data (42 subjects x 19 conditions cell array (246 chans x 4 freq bands))
% flatten the data in each cell

flatten_by_band = @(x) x(:).'; % preserve channels' psds for each band 
% (transpose for row instead of column representation)

function flattened_data = flat_by_channel(x) % preserve band values for each channel
    x_t = transpose(x);
    flattened_data = x_t(:).'; 
end

flatten_by_channel = @(x) flat_by_channel(x); %%%%%%%%% try this

data = cellfun(flatten_by_channel, data, "UniformOutput", false); %flatten_by_band 

trials_features = [];
% concatenate the flattened data in a vertical manner in an evoked
% trial/condition per subject * feauture manner
for r=1:size(data,1)
    for c=1:size(data,2)
        trials_features = [trials_features; data{r,c}];
    end
end

% concatenate the unique labels for the 42 subjects in a column
labels = [];
for sub=1:size(data, 1)
    labels = [labels u_labels];
end 

%% set labels for classification: %%

% remove unused classes
function [labels, trial_matrix] = remove_class(trial_matrix, orig_labels, classes_toRmv)
    trial_matrix(ismember(orig_labels, classes_toRmv), :) = [];
    labels = removecats(orig_labels, classes_toRmv); % sets classes_toRmv --> <undefined> but doesn't delete the space from the array
    labels = labels(ismember(labels, categories(labels)));
end

tic 

classes_toRmv = {'odd'};
[labels, trials_features] = remove_class(trials_features, labels, classes_toRmv);

toc

%% F1-F2 classification
mergeLabels = mergecats(labels, {'FS1', 'FM1', 'FL1'}, 'F1');
mergeLabels = mergecats(mergeLabels, {'FS2', 'FM2', 'FL2'}, 'F2');

%% 2 class classifcation by repetition

mergeLabels = mergecats(labels, {'FS1', 'FM1', 'FL1', 'PS1', 'PM1', 'PL1', 'NS1', 'NM1', 'NL1'}, 'Rep_1');
mergeLabels = mergecats(mergeLabels, {'FS2', 'FM2', 'FL2', 'PS2', 'PM2', 'PL2', 'NS2', 'NM2', 'NL2'}, 'Rep_2');

%% 3 class classification by category

mergeLabels = mergecats(labels, {'FS1', 'FM1', 'FL1', 'FS2', 'FM2', 'FL2'}, 'Food');
mergeLabels = mergecats(mergeLabels, {'PS1', 'PM1', 'PL1', 'PS2', 'PM2', 'PL2'}, 'Positive');
mergeLabels = mergecats(mergeLabels, {'NS1', 'NM1', 'NL1', 'NS2', 'NM2', 'NL2'}, 'Neutral');

%% 3 class classification by delay
mergeLabels = mergecats(labels, {'FS1', 'FS2', 'PS1', 'PS2', 'NS1', 'NS2'}, 'Short');
mergeLabels = mergecats(mergeLabels, {'FM1', 'FM2', 'PM1', 'PM2', 'NM1', 'NM2'}, 'Medium');
mergeLabels = mergecats(mergeLabels, {'FL1', 'FL2', 'PL1', 'PL2', 'NL1', 'NL2'}, 'Long');
%% 6 class classification by category*repetition

mergeLabels = mergecats(labels, {'FS1', 'FM1', 'FL1'}, 'F1');
mergeLabels = mergecats(mergeLabels, {'FS2', 'FM2', 'FL2'}, 'F2');
mergeLabels = mergecats(mergeLabels, {'PS1', 'PM1', 'PL1'}, 'P1');
mergeLabels = mergecats(mergeLabels, {'PS2', 'PM2', 'PL2'}, 'P2');
mergeLabels = mergecats(mergeLabels, {'NS1', 'NM1', 'NL1'}, 'N1');
mergeLabels = mergecats(mergeLabels, {'NS2', 'NM2', 'NL2'}, 'N2');

%% 6 class classification by delay*repetition

mergeLabels = mergecats(labels, {'FS1', 'PS1', 'NS1'}, 'S1');
mergeLabels = mergecats(mergeLabels, {'FS2', 'PS2', 'NS2'}, 'S2');
mergeLabels = mergecats(mergeLabels, {'FM1', 'PM1', 'NM1'}, 'M1');
mergeLabels = mergecats(mergeLabels, {'FM2', 'PM2', 'NM2'}, 'M2');
mergeLabels = mergecats(mergeLabels, {'FL1', 'PL1', 'NL2'}, 'L1');
mergeLabels = mergecats(mergeLabels, {'FL2', 'PL2', 'NL2'}, 'L2');

%% 9 class classification by category*delay

mergeLabels = mergecats(labels, {'FS1','FS2'}, 'FS');
mergeLabels = mergecats(mergeLabels, {'FM1','FM2'}, 'FM');
mergeLabels = mergecats(mergeLabels, {'FL1','FL2'}, 'FL');
mergeLabels = mergecats(mergeLabels, {'PS1','PS2'}, 'PS');
mergeLabels = mergecats(mergeLabels, {'PM1','PM2'}, 'PM');
mergeLabels = mergecats(mergeLabels, {'PL1','PL2'}, 'PL');
mergeLabels = mergecats(mergeLabels, {'NS1','NS2'}, 'NS');
mergeLabels = mergecats(mergeLabels, {'NM1','NM2'}, 'NM');
mergeLabels = mergecats(mergeLabels, {'NL1','NL2'}, 'NL');
%% 18 class classification

mergeLabels = labels;

%% extract number of unique classes 

n_classes = length(categories(mergeLabels));

%% split data to training and testing sets
% seed for cvpartition, for reproducibility
p=0.2; % fraction of test data
training_test_split = cvpartition(mergeLabels, "HoldOut", p);
train_set = trials_features(training_test_split.training, :);
train_labels = mergeLabels(training_test_split.training);
test_set = trials_features(training_test_split.test, :);
test_labels = mergeLabels(training_test_split.test); 

%% PCA
[coeff,train_transformed,~,~,explained,mu] = pca(train_set);
test_transformed = (test_set-mu)*coeff;

%% Simple model with cross validation 
c = cvpartition(train_labels, KFold=5, Stratify=true);
trained_model = fitcecoc(train_transformed, train_labels, 'Learners', 'svm', 'OptimizeHyperparameters', ...
    'all', 'HyperparameterOptimizationOptions', struct('Useparallel', true, 'CVPartition', c, 'MaxObjectiveEvaluations', 100));

%% Prediction and confusion matrix on test data
% Predict the labels of the test data using the trained SVM model
[predicted_labels, scores] = predict(trained_model, test_transformed);

% create a confusion matrix based on true and predicted labels
% 'row-normalized' - Normalize each cell value by the number of observations that has the same true class.
figure()
confusionMat = confusionchart(test_labels, predicted_labels, 'Normalization', 'row-normalized', 'Title', 'Confusion matrix of decoder preformance');

accuracy = sum(predicted_labels.' == test_labels)/numel(test_labels)*100;
fprintf('Classifier accuracy: %.2f', accuracy)

clearvars("trained_model")
