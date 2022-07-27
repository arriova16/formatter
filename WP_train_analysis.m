%%% New formatting script
data_folder = 'Z:\Tori\ICMS\data\7_25-7_29\WP'
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
        if ismissing(temp_data{2,2})
            % Remove 2nd column
            temp_data = temp_data(:,[1,3:end]);
        elseif ismissing(temp_data{2,15})
            temp_data(1,[2:14]) = temp_data(1,[3:15]);
            temp_data = temp_data(:,1:end-1);
        end
        
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
    block_struct(ii).TestElec = response_table{1,13};
    block_struct(ii).ResponseTable = response_table;
    block_struct(ii).NumTrials = size(response_table,1);    

    ii = ii + 1;
end

%% WP percent correct

%for t = 1:length(block_struct)

%correct_total = rdivide(250,block_struct(t).NumTrials);

%percent_correct = correct_total * 100;

%block_struct(ii).Percent_Correct = percent_correct;

%end

