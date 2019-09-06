clear all;

data_type = '%g';
data_separator = ',';
dist_type = 'euclidean';
ignore_error = 1; % this is used in distance calculation(affinity) where few error bits of very high value will be ignored.

%% Get testing data i.e. Antigens
%<<
    [data file_used] = get_data(data_type, data_separator, 'Testing File','test_*.data');
    
    
    
    %<< process raw data
% % % %     total_features = size(data,2) - 1;
% % % %     file_format = '';
% % % %     for(i = 1:total_features)
% % % %         file_format = [file_format data_type data_separator];
% % % %     end
% % % %     file_format = [file_format data_type];
% % % %     file_format = [file_format '\n'];
% % % % 
% % % %     norm_data = norm_variance(data(:,[1:total_features])); 
% % % %     norm_data = norm_scale01(norm_data);    
% % % %     data = [norm_data data(:,total_features+1)];
% % % % 
% % % %     min_val = 0; % Since scaled to [0 1]
% % % %     max_val = 1;
    %>>
    
    
    
    
    Ag = data; % testing data antigens
    M = size(Ag,1); %Antigen population size
    test_classifications = Ag(:,size(Ag,2));
%>>

%% Get memory cells 
%<<
    [data file_used] = get_data(data_type, data_separator, 'Memory File','*data_memory.out');
    Mem = data; % memory cells + class
    m = size(Mem,1); % size of memory cells
    L = size(Mem,2) - 1; %Antigen's epitope's length => dimensions/features
%>>

Ab_struct = clonal_struct('Antibody');
Ag_struct = clonal_struct('Antigen');
            
for(i = 1:M)
    Ag_struct(i).type = 'antigen';
    Ag_struct(i).epitope = Ag(i,[1:L]); % L are features - Antigen file can be labelled or non labelled
    Ag_struct(i).class = -1; % objective is to do the classification
end

clear Ag;            
Ag = Ag_struct;
clear Ag_struct;

for(i = 1:m)  
    Ab_struct(i).type = 'anitibody';
    Ab_struct(i).receptor = Mem(i,[1:L]); 
    Ab_struct(i).class = Mem(i,L+1); % The last column has class number
    Ab_struct(i).affinity = -1; %>>Mem(i,L+2); % NOT NEEDED
    Ab_struct(i).isMemoryCell = true;
end

Ab = Ab_struct;
clear Ab_struct;

temp1 = reshape([Ag(:).epitope],[L M])';
temp2 = reshape([Ab(:).receptor],[L m])';
temp = [temp1' temp2']';
max_dist = maximum_possible_distance( temp ); % maximum distance in given vicinity

clear temp1 temp2 temp;

classification_threshold = 0.1; % 79%

for(i = 1:M)
    f(i,:) = affinity(reshape([Ab(:).receptor], [L m])', [Ag(i).epitope], dist_type , max_dist, ignore_error);
    [aff position] = sort(f(i,:),'descend'); % affinity ranks of antibodies

    if(aff(1) > aff(2)) % only one best value, 
        %%make sure both belongs to diffent class

        if(aff(1) > classification_threshold)
            Ag(i).class = Ab(position(1)).class;
        else
            Ag(i).class = -1; % unable to classify
        end
    else %more than one same best val
        best_vals_idx = find(aff == aff(1)); % first k elements            
        total_best_vals = size(best_vals_idx,1);

        % if all same best vals belongs to same class then easilty classify
        % to the same class.
        all_best = [Ab(position([best_vals_idx])).class];
        is_all_best = find(all_best == Ab(position(best_vals_idx(1))).class);

        if(length(is_all_best) == length(all_best))
            if(aff(1) > classification_threshold)
                Ag(i).class = Ab(position(1)).class;
            else
                Ag(i).class = -1; % unable to classify
            end
        else
            %find(aff > classification_threshold);
            Ag(i).class = -2; % overlap - classification. suggestion: increase threshold    
        end                        
    end
end

fprintf('Tesing Completed!\n\n');

out = [[Ag(:).class]' test_classifications];

fprintf('Preidicted Class ---- Given Class\n');
fprintf('      %d                    %d\n',out');

TP = 0;
for(i = 1: M)
    if(Ag(i).class == test_classifications(i))
        TP = TP+1;
    end
end

fprintf('True Poisitive %1.2f\n',TP/M*100);
