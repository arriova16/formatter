%%% New formatting script
%data_folder = '\\BENSMAIA-LAB\LabSharing\Tori\Testing\Test';
%data_folder =  '\\BENSMAIA-LAB\LabSharing\Tori\ICMS\data\5_31-6_3\pinot\training';
data_folder = 'Z:\Tori\ICMS\data\7_25-7_29'
file_list = dir(data_folder);

%% Load and parse .rsp files
ii = 1;
block_struct = struct();
for i = 1:length(file_list)
    % Check if file is .rsp
    if ~contains(file_list(i).name, '.rsp') || ~contains(file_list(i).name, 'discrimination')
        continue
    end
    
    % Parse filename
    fname_split = strsplit(file_list(i).name, '_');
    %data_struct(ii).StudyID = fname_split{1};
    block_struct(ii).Animal = fname_split{2}(1:end-7);
    block_struct(ii).Protocol = fname_split{3};
    t_idx = find('20220223T081143' == 'T');
    block_struct(ii).Date = datestr(datenum(fname_split{4}(1:t_idx-1), 'yyyymmdd'), 1);
    %data_struct(ii).Time = datestr(datenum(fname_split{4}(t_idx+1:end-4), 'hhmmss'), 13);
    
    % Load the data
    %temp_data = readmatrix(fullfile(data_folder, file_list(i).name), 'Delimiter', 't','FileType','text', 'NumHeaderLines', 1);
    temp_data = readcell(fullfile(data_folder, file_list(i).name), 'FileType', 'text');
    if size(temp_data,2) == 15
        % Remove 2nd column
        temp_data = temp_data(:,[1,3:end]);
        
    elseif size(temp_data,2) == 24
        temp_data(1,2:end-3) = temp_data(1,[3:11,13:20,22:end]);
        temp_data = temp_data(:,[1:4,8:12,16:20]);
        cn = [4:13];
        rep_string = {'CondStimDur','CondStimAmp','CondStimFreq','CondStimPhase',...
                      'CondElectrodeNum','TestStimDur','TestStimAmp','TestStimFreq',...
                      'TestStimPhase','TestElectrodeNum'};
        for j = 1:length(cn)
            temp_data{1,cn(j)} = rep_string{j};
        end
    end
    % Remove aborted trials
    abort_idx = strcmpi(temp_data(:,14), 'empty') | strcmpi(temp_data(:,14), 'no');
    temp_data = temp_data(~abort_idx,:);    
    % Remove footer
    foot_idx = cellfun(@(c) isnumeric(c), temp_data(2:end,1));
    temp_data = temp_data([true;foot_idx],:);
    % Convert to table
    response_table = cell2table(temp_data(2:end,:), 'VariableNames', temp_data(1,:));
    block_struct(ii).CondElec = response_table{1,8};
    block_struct(ii).TestElec = response_table{1,13};
    if block_struct(ii).CondElec~= block_struct(ii).TestElec
        block_struct(ii).Period = 'Pre';
    else
        block_struct(ii).Period = 'Post';
    end
    block_struct(ii).ResponseTable = response_table;
    block_struct(ii).NumTrials = size(response_table,1);    

    ii = ii + 1;
end

%% Get p(detected) & d' for each block
for i = 1:length(block_struct)
    [u_test_amps, ~, ia] = unique(block_struct(i).ResponseTable.TestStimAmp);
    p_detect = zeros([length(u_test_amps),1]);
    for j = 1:length(u_test_amps)
        correct_idx = strcmp(block_struct(i).ResponseTable.Response(ia == j), 'correct');
        p_detect(j) = sum(correct_idx) / length(correct_idx);
        
    end
    % Correct for catch trials
    if any(u_test_amps == 0)
        catch_idx = find(u_test_amps == 0);
        p_detect(catch_idx) = 1 - p_detect(catch_idx);
    end
    
    % Compute d' from p(detect) where  d' = norminv(p(hit)) - norminv(p(false alarm))
    dprime = NaN([length(u_test_amps),1]);
    pmiss = p_detect(1);
    if pmiss == 0 % Correct for 0 false alarm
        pmiss = 0.001;
    end
    for j = 1:length(dprime)-1
        phit = p_detect(j+1);
        if phit == 1 % Correct for infinite hit rate
            phit = .999;
        end
        dprime(j+1) = norminv(phit) - norminv(pmiss);        
    end
    
    % Make a table & add to struct
    block_struct(i).DetectionRates = array2table([u_test_amps, p_detect, dprime], 'VariableNames', {'Amplitude', 'pDetect', 'dPrime'});
end


%% Day by day psychometric plots
sigfun = @(c,x) 1./(1 + exp(-c(1).*(x-c(2)))); % c(1) = rate of change, c(2) = x-offset
opts = optimset('Display','off');
xq = linspace(0,60);

u_animals = unique({block_struct.Animal});
%u_dates = unique([{data_struct.Date}', '], 'rows');
for a = 1:length(u_animals)
    a_idx = strcmp({block_struct.Animal}, u_animals{a});
    animal_dates = unique({block_struct(a_idx).Date});
    for i = 1:length(animal_dates)
        figure('Name', sprintf('%s - %s - Ch%d',u_animals{a}, animal_dates{i}, block_struct(a_idx).TestElec));
        % Find the date & pre
        pre_idx = find(a_idx &...
                   strcmp({block_struct.Date}, animal_dates{i}) &...
                   strcmp({block_struct.Period}, 'Pre'));
        if ~isempty(pre_idx)
            sp1 = subplot(1, 2, 1, 'parent', gcf); hold on
            scatter(block_struct(pre_idx).DetectionRates{:,1}, block_struct(pre_idx).DetectionRates{:,2},...
                'MarkerFaceColor', [.6 .6 .6],'MarkerEdgeColor', [.6 .6 .6], 'MarkerFaceAlpha', .2, 'Parent', sp1)
            [sig_c, resnorm] = lsqcurvefit(sigfun,[1.7/30,30], block_struct(pre_idx).DetectionRates{:,1},...
                                                               block_struct(pre_idx).DetectionRates{:,2},...
                                                               [0, 0],[1, 60], opts);
            plot(xq, sigfun(sig_c, xq), 'Color', [.6 .6 .6], 'Parent', sp1)
            xticks([0:10:60]); 
    
            sp2 = subplot(1, 2, 2); hold on
            if size(block_struct(pre_idx).DetectionRates,1) < 3
                scatter(block_struct(pre_idx).DetectionRates{2,1}, block_struct(pre_idx).DetectionRates{2,3},...
                'MarkerFaceColor', [.6 .6 .6],'MarkerEdgeColor', [.6 .6 .6], 'MarkerFaceAlpha', .2, 'Parent', sp2)
            else
                plot(block_struct(pre_idx).DetectionRates{2:end,1}, block_struct(pre_idx).DetectionRates{2:end,3}, 'Color', [.6 .6 .6],...
                    'Parent', sp2)
            end
        end

        % Find the date & post
        post_idx = find(a_idx &...
                   strcmp({block_struct.Date}, animal_dates{i}) &...
                   strcmp({block_struct.Period}, 'Post'));
        if ~isempty(post_idx)
            scatter(block_struct(post_idx).DetectionRates{:,1}, block_struct(post_idx).DetectionRates{:,2},...
                'MarkerFaceColor', [1 .6 .6],'MarkerEdgeColor', [1 .6 .6], 'MarkerFaceAlpha', .2, 'Parent', sp1)
            [sig_c, resnorm] = lsqcurvefit(sigfun,[1.7/30,30], block_struct(post_idx).DetectionRates{:,1},...
                                                               block_struct(post_idx).DetectionRates{:,2},...
                                                               [0, 0],[1, 60], opts);
            plot(xq, sigfun(sig_c, xq), 'Color', [1 .6 .6], 'Parent', sp1);
            xticks([0:10:60]);

            if size(block_struct(post_idx).DetectionRates,1) < 3
                scatter(block_struct(post_idx).DetectionRates{2,1}, block_struct(post_idx).DetectionRates{2,3},...
                'MarkerFaceColor', [1 .6 .6],'MarkerEdgeColor', [1 .6 .6], 'MarkerFaceAlpha', .2, 'Parent', sp2)
            else
                plot(block_struct(post_idx).DetectionRates{2:end,1}, block_struct(post_idx).DetectionRates{2:end,3}, 'Color', [1 .6 .6],...
                    'Parent', sp2)
            end
        end

        ylabel(sp1, 'p(Detected)')
        ylabel(sp2, 'd''')

        xlabel(sp1,sprintf('Amplitude (%sA)', GetUnicodeChar('mu')))
        xlabel(sp2,sprintf('Amplitude (%sA)', GetUnicodeChar('mu')))
    end
end

%% Combine across days
channel_struct = struct();
unique_electrodes = struct2table(block_struct);
unique_electrodes = unique_electrodes(:, ["Animal", "TestElec"]);
unique_electrodes = unique(unique_electrodes, 'rows'); % Remove period duplicates

% Find index matching unique values eg line 112
for i = 1:size(unique_electrodes,1)
    % Assign unique to new struct
    channel_struct(i).Animal = unique_electrodes.Animal{i};
    channel_struct(i).Channel = unique_electrodes.TestElec(i);

    % Find matching from block struct
    pre_idx = find(strcmp({block_struct.Animal}, unique_electrodes.Animal{i}) &...
               [block_struct.TestElec] == unique_electrodes.TestElec(i) &...
               strcmp({block_struct.Period}, 'Pre'));

    post_idx = find(strcmp({block_struct.Animal}, unique_electrodes.Animal{i}) &...
               [block_struct.TestElec] == unique_electrodes.TestElec(i) &...
               strcmp({block_struct.Period}, 'Post'));
    
    % Combine pre-sessions
    if ~isempty(pre_idx)
        channel_struct(i).ControlRT = cat(1, block_struct(pre_idx).ResponseTable);
    end
    
    % Combine post-sessions
    if ~isempty(post_idx)
        channel_struct(i).BoostRT = cat(1, block_struct(post_idx).ResponseTable);
    end
end

%% dprime and pdetect of controlrt and boost rt
for i = length(channel_struct)
    [u_test_amps, ~, ia] = unique(channel_struct(i).BoostRT.TestStimAmp);
    p_detect_comb = zeros([length(u_test_amps),1]);
    for j=1:length(u_test_amps)
        correct_idx_comb = strcmpi(channel_struct(i).BoostRT.TestStimAmp);
        p_detect_comb(j) = sum(correct_idx_comb)/ length(correct_idx_comb);
    end
    if any(u_test_amps == 0)
        catch_idx_comb = find(u_test_amps_comb == 0);
        p_detect_comb(catch_idx_comb) = 1 - p_detect_comb(catch_idx_comb);
    end


end


%% Get p(detected) & d' for each block
%for i = 1:length(block_struct)
    %[u_test_amps, ~, ia] = unique(block_struct(i).ResponseTable.TestStimAmp);
    %p_detect = zeros([length(u_test_amps),1]);
    %for j = 1:length(u_test_amps)
        %correct_idx = strcmp(block_struct(i).ResponseTable.Response(ia == j), 'correct');
        %p_detect(j) = sum(correct_idx) / length(correct_idx);
        
    %end
    %% Correct for catch trials
    %if any(u_test_amps == 0)
        %catch_idx = find(u_test_amps == 0);
        %p_detect(catch_idx) = 1 - p_detect(catch_idx);
    %end
    
    % Compute d' from p(detect) where  d' = norminv(p(hit)) - norminv(p(false alarm))
    %dprime = NaN([length(u_test_amps),1]);
    %pmiss = p_detect(1);
    %if pmiss == 0 % Correct for 0 false alarm
        %pmiss = 0.001;
    %end
   %for j = 1:length(dprime)-1
        %phit = p_detect(j+1);
        %if phit == 1 % Correct for infinite hit rate
            %phit = .999;
        %end
        %dprime(j+1) = norminv(phit) - norminv(pmiss);        
    %end
    
    % Make a table & add to struct
    %block_struct(i).DetectionRates = array2table([u_test_amps, p_detect, dprime], 'VariableNames', {'Amplitude', 'pDetect', 'dPrime'});
%end


