% Developed by Anurag Sharma

% version 3.5 - introduced noise removal code

% version 3.4 - I am thinking of putting penalty-reward feature for
% classification. If attraction towards opposite class is high then drop
% that memory cell.

% NOTE: THIS ALGORITHM ALWAYS "MUST" PREPROCESS THE DATA TO SCALE IT TO    
% [0 1] RANGE. AND THIS RANGE IS USED IN MANY SUBROUTINES....

% Version 3.3 - those antignes which are considered will be removed from
% repertoire for that generation.
% Note: duplication checking in memory repertoire has been removed as it
% slows down the process.
% POSSIBLE DANGER (allow duplication): all memeory cells may converge to
% one cell... becuase it removes the min from repertoire.

% Version 3.2 - consider many best clones per Antigen
% Version 3.1 - Updated version 3.0 - solving memory usage problem
% Version 3.0 - very expensive... problem with memory usage and speed.
% This one has m < M
% IDEA: One memory cell can be combined with more than one antigen.
% It will create generalized memory cells
% since m < M hence m is divided into few slots. each slot belongs to one
% class. sizes for each slot is of the same ratio as the size of classes in
% training antigen of size M.
% only best clone for a given antigen is picked
% then its average affinity with all other antigens is calculated to make a
% good generalized memory cell
% if this average affinity is greater than the minimum affinity in the
% class slot then it will conditionally replace the minimum affinity. It
% happens if there is no exact copy of memory already in the slot.

clc;
clear all;

% % % %% data preprocessing
% % % %<<
    %step 1: describe data
    data_type = '%g';
    data_separator = ',';
    dist_type = 'euclidean';
    mutation_type = 'euclidean';
    
%         data_type = '%d';
%     data_separator = ' ';
%     dist_type = 'hamming';
%     mutation_type = 'bit flip';
    
    
    %step 2: get raw data
    [data file_used] = get_data(data_type,data_separator,'Raw Data File','*.data');
    data_size = size(data,1);
    total_features = size(data,2) - 1;
    
    %step 3: normalize data
    %can also use norm_variance(data); % to make variance = 1
% % %     norm_data = norm_variance(data(:,[1:total_features])); 
% % %     norm_data = norm_scale01(norm_data);
    
    for (i = [1:total_features])
        norm_data(:,i) = norm_variance(data(:,i));
        norm_data(:,i) = norm_scale01(norm_data(:,i));
    end
    
    
    data = [norm_data data(:,total_features+1)];
    clear norm_data;
    min_val = 0; % Since scaled to [0 1]
    max_val = 1; % Since scaled to [0 1]
    
    %step 4: separate training and testing data
    training_ratio = 0.8;
    
    sorted_data = sortrows(data, total_features + 1); % sort according to class
    class_list = unique(data(:,total_features + 1));
    
    count = 0;
    for (i = [class_list]')
        count = count + 1;
        class_size(count) = length(find (data(:,total_features + 1) == i));        
    end    
    
    cur_pointer = 0;
    train_data = [];
    test_data = [];
    
    for (i = 1:count)
        rand_order = cur_pointer + randperm(class_size(i));
        cur_pointer = cur_pointer + class_size(i);
        
        training_size = round(class_size(i) * training_ratio);
        testing_size = class_size(i) - training_size;
        train_data = [train_data' sorted_data(rand_order(1:training_size),:)']';
        test_data = [test_data' sorted_data(rand_order(training_size+1:class_size(i)),:)']';    
    end

    %<<    
%     training_size = round(data_size * training_ratio);
%     testing_size = data_size - training_size;
%     train_data = [data([1:training_size],:)']';
%     test_data = [data([training_size+1:data_size],:)']';    
% % %     for (i = 1:count)
% % %         rand_order = cur_pointer + [1:class_size(i)];
% % %         cur_pointer = cur_pointer + class_size(i);
% % %         
% % %         training_size = round(class_size(i) * training_ratio);
% % %         testing_size = class_size(i) - training_size;
% % %         train_data = [train_data' sorted_data(rand_order(1:training_size),:)']';
% % %         test_data = [test_data' sorted_data(rand_order(training_size+1:class_size(i)),:)']';   
% % %     end
    %>>
    
    %step 5: print training and testing data into files
    train_file = ['data\train_' file_used];
    test_file = ['data\test_' file_used];
    
    file_format = '';
    for(i = 1:total_features)
        file_format = [file_format data_type data_separator];
    end
    file_format = [file_format data_type];
    file_format = [file_format '\n'];
    
    train_out = fopen(train_file,'w');    
    fprintf(train_out, '%d \n\n', total_features + 1);
    fprintf(train_out,file_format, train_data');
    fclose(train_out);
    
    test_out = fopen(test_file,'w');    
    fprintf(test_out, '%d \n\n', total_features + 1);
    fprintf(test_out,file_format, test_data');
    fclose(test_out);
%>>
  
clear data;
data = train_data;
clear train_data test_data;

M = size(data,1); %Antigen population size
L = size(data,2) - 1; %Antigen's epitope's length => dimensions/features

% Randomly generate antibodies
Ngen = 8; % was 70;
%m = M; % older version: it is assumed m >= M ... problem ... it is always m=M
m = 100; % max(30,round(M*0.5)); % current version(2): m < M 
r = round(0.10 * m); % 10% of m
d = 0; %round(0.50 * r); % d < r 50% of r
N = m + r; % antibody population
ro = 5;
n = 30; %max(10, round(m*0.15)); % n best Antibodies(with higher affinities) to be selected.
k = 10; %min(n, round(m*0.05)); %n; % or simply make k = n; k best clones to be selected
min_memory_affinity = 0.7;

%%%% IN case if none antibodies have "good enough" affinity this algorithm
%%%% is still forced to select n best.... *****************  THIS SECTION
%%%% CAN BE IMROVED................
beta = 0.5; % multiplying factor 
init_train_memory_ratio = 8; %floor(M/m); % maximum number of antigens per memory cell.
% min_neighborhood_size = 3;

ignore_error = 0.9; % this is used in distance calculation(affinity) where few error bits of very high value will be ignored.

Ab_struct = clonal_struct('Antibody');        
Ag_struct = clonal_struct('Antigen'); 
            
for(i = 1:M)
    Ag_struct(i).type = 'antigen';
    Ag_struct(i).epitope = data(i,[1:L]); % 1-L are features
    Ag_struct(i).class = data(i,L+1); %The last column L+1 denotes class number
end
       
Ag = Ag_struct;
clear Ag_struct;
Ab = generate_random_numbers(N, L, min_val, max_val, data_type);

% temp1 = reshape([Ag(:).epitope],[L M])';
% temp2 = Ab;
% temp = [temp1' temp2']';
% max_dist = maximum_possible_distance(temp); % maximum distance in given vicinity
% clear temp1 temp2 temp;

UpperBound = 0;
for(i = 1:L) % since all are normalized to [0, 1]
    UpperBound = UpperBound + 1*1;
end
UpperBound = sqrt(UpperBound); % NOTE: this is not same as max_dist it is more than that

for(i = 1:N)
    Ab_struct(i).type = 'anitibody';
    Ab_struct(i).receptor = Ab(i,:); 
    Ab_struct(i).class = -1; % intially not assigned to any class
    Ab_struct(i).affinity = -1; % intially no affinity 
    Ab_struct(i).isMemoryCell = false;
end

for(i = r+1:N)
    Ab_struct(i).isMemoryCell = true;
end

clear Ab;
Ab = Ab_struct;
clear Ab_struct;

class_info = clonal_struct('class_info');
[class_info Ab] = update_class_data(Ag, class_info, Ab, m, r, M);

% Ab(1:N,L+1) = 0; % second last column is indicator whether it is a memory cell
% Ab(1:N,L+2) = -1; % last column is an indicator for which anitigen it is bounded with
% assume first {r} elements are for set Ab{r} and second {m} elements are
% for set Ab{m}.

class_size = size(class_info,2);

% % % last_cell = r;
% % % for ( i = 1:class_size)
% % %     rand_Ab = last_cell + randperm(class_info(i).size);
% % %     
% % %      
% % %     this_Ab_class_size = size(class_info(i).Antibody_NRange,2);
% % %     
% % %     counter = 1;
% % %     for (temp = [class_info(i).Antibody_NRange])
% % %         Ab(temp).receptor = Ag(rand_Ab(counter)).epitope;
% % %         counter = counter + 1;
% % %     end
% % %     
% % %     last_cell = max(class_info(i).locations);
% % % end




tic;
tend = 0;

init_Ag = Ag;
init_M = M;
init_class_info = class_info;



for t = 1: Ngen
    Ag = init_Ag;
    M = init_M;
    class_info = init_class_info;
    j = 1;
    while(j <= M)    % step 1        
        
        % step 2: get affinity for all types of antibodies. m + r
        Ag_class_j = Ag(j);
        
        %<< temporary measure
% % %             if(t > 8 && Ag_class_j.class == 1)
% % %                 j = j + 1;
% % %                 
% % %                 
% % %                 init_train_memory_ratio = [5 2];
% % %                 continue;
% % %             end
        %>>
        
        const_class_info_idx = find([init_class_info(:).number] == Ag_class_j.class);
        const_Ag_j_class_mapped_NRange = init_class_info(const_class_info_idx ).Antibody_NRange;
        Abs_used = [[1:r] const_Ag_j_class_mapped_NRange];
        size_Abs_used = size(Abs_used,2);
        train_memory_ratio = init_train_memory_ratio;
        
        f_j = affinity(reshape([Ab(Abs_used).receptor], [L size_Abs_used])', [Ag_class_j.epitope], dist_type, UpperBound);
        [aff position] = sort(f_j,'descend');
        
        % step 3
% % %         aff_used = aff(1:n_max);
% % %         far_away_abs = find ( [aff_used] < (mean(aff_used) - std(aff_used)));          
% % %         aff_used(far_away_abs) = [];
% % %         n = length(aff_used);
        Abn = reshape([Ab([ Abs_used(position(1:n))]).receptor],[L n])';  %select best n antibodies(i.e. antibodies with n highest affinity values)                    
        %%%%% IDEA: one anitbody can be good for one antigen but bad for
        %%%%% other.... we should not remove an antibody based on just one
        %%%%% particular antigen........... think......
        
        
        %step 4        
        [C_j{1,1} clone_sizes] = clone(Abn, beta, size_Abs_used); % cloning        
% %         clone_sizes
        
        % step 5
        C_mutated_j{1,1} = hypermutate(aff(1:n),ro,C_j{1,1},L, mutation_type, min_val, max_val); %C{j,1}{:,1}
        clear C_j;
        
        % step 6
        Ab_mutated = cell2mat(C_mutated_j{1,1}); 
        clear C_mutated_j;
        f_mutated_j = affinity(Ab_mutated, [Ag_class_j.epitope], dist_type, UpperBound);
        [best_affs best_positions] = sort(f_mutated_j,'descend');                
        
        %[best_affs uniq_locs] = unique(best_affs); % in ascending order i.e. lower affinity to higher        
        %size_unique_bests = size(best_affs,2);        
        %best_k = min(k,size_unique_bests); 
        
        best_k = k; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        if(best_k > length(f_mutated_j))
            best_k = length(f_mutated_j);
        end
        
        
        class_info_idx = find([class_info(:).number] == Ag_class_j.class);
        Ag_j_class_locations = [class_info(class_info_idx ).locations];
        Ag_j_class_size = class_info(class_info_idx ).size;
        
        k_counter = 0;
        while(k_counter < best_k) 
            k_counter = k_counter+1;
            %best_aff = best_affs(size_unique_bests - i + 1);  % get the first one, second one... up to best ith one 
            %best_position = best_positions(uniq_locs(size_unique_bests - i + 1)); %get only the ith best, NOT n antibodies
            best_aff = best_affs(k_counter);  % get the first one, second one... up to best ith one 
            best_position = best_positions(k_counter); %get only the ith best, NOT n antibodies
            
            Ab_mutated_best = Ab_mutated(best_position,:);  %select the best antibody(highest affinity)
            % we will keep the best antibodies for each antigens
            % at the end we will have M best antibodies for M antigens. one
            % antibody for one antigen.


            
            %ANU step(a): Get average affinity
            % Compare best mutated cell with all "antigens of same class"
            % and calculate avg.        

            
            if(best_k > Ag_j_class_size)
                best_k = Ag_j_class_size;
            end
        
            Ag_j_class_f_j = affinity(reshape([Ag(Ag_j_class_locations).epitope], [L Ag_j_class_size])', Ab_mutated_best, dist_type, UpperBound);
            [Ag_j_class_aff Ag_j_class_position] = sort(Ag_j_class_f_j,'descend');        
            % get top train_memory_ratio        
            
            if (train_memory_ratio > min(size(Ag_j_class_aff,2)))
                error('error!! train_memory_ratio > min(size(Ag_j_class_aff,2))');
                break;
            end
            %%%train_memory_ratio = min(size(Ag_j_class_aff,2),train_memory_ratio);            
            
            best_Ag_j_class = Ag_j_class_aff(1:train_memory_ratio);                      

%             [tval tloc] = max(abs(diff([best_Ag_j_class])));
%             min_aff_considered = best_Ag_j_class(tloc);
        

            
% % %             far_neighbors = find ( [best_Ag_j_class] < (mean([best_Ag_j_class]) - std([best_Ag_j_class])));          
% % %             best_Ag_j_class(far_neighbors) = [];
            
            size_best_Ag_j_class = length(best_Ag_j_class);
            
            avg_aff = mean([best_Ag_j_class]);  % average affinity of best mutated memory cell   
            min_aff = min([best_Ag_j_class]);
            stdev = std([best_Ag_j_class]);
            
%             if( lower_bound < train_memory_ratio )
%                 lower_bound
%             end
            
%             if (tloc <= 3)
%                 continue;
%             end
               

            Ag_j_class_mapped_NRange = class_info(class_info_idx ).Antibody_NRange;
            cur_aff = [Ab(Ag_j_class_mapped_NRange).affinity]; %% IMP I don't know why it slows the process.????????????????????????????????????????????

    %         Ab(Ag_j_class_mapped_NRange).affinity = ? I DON'T THINK THAT IT
    %         IS NEEDED

            [min_class_aff min_class_loc] = min(cur_aff);
            min_class_loc_in_NRange = Ag_j_class_mapped_NRange(1) - 1 + min_class_loc;

            %check if this mutated cell is not already in Ab
            %<<
% %             mem_class_size = size(Ag_j_class_mapped_NRange,2);
% %             mem_set = reshape([Ab(Ag_j_class_mapped_NRange).receptor], [L mem_class_size])';
% %             str2 = num2str(Ab_mutated_best);       
% % 
            bexist = false;
% %             for(mem_cell = 1:mem_class_size)
% %                 str1 = num2str(mem_set(mem_cell,:));
% %                 if(strcmp(str1(1,:) , str2(1,:)))
% %                     bexist = true;
% %                     break;
% %                 end
% %             end
            %>>

            
             % Check average affinity with other classes....
            % If this attracts more towards other class(es) than reject it.
            % can be moved up.....
            
            if(~bexist)
                OTHER_class_info_idx = find([class_info(:).number] ~= Ag_class_j.class);            
                for(o_c = [OTHER_class_info_idx])
                    OTHER_Ag_j_class_locations = [class_info(o_c ).locations];
                    OTHER_Ag_j_class_size = class_info(o_c ).size;
                    OTHER_Ag_j_class_f_j = affinity(reshape([Ag(OTHER_Ag_j_class_locations).epitope], [L OTHER_Ag_j_class_size])', Ab_mutated_best, dist_type, UpperBound);

% %                     [OTHER_Ag_j_class_aff OTHER_Ag_j_class_position] = sort(OTHER_Ag_j_class_f_j,'descend');  % change this to get max... will be faster... 
                    OTHER_Ag_j_class_aff = max(OTHER_Ag_j_class_f_j);
% %                     OTHER_train_memory_ratio = min(size(OTHER_Ag_j_class_aff,2),train_memory_ratio);
% %                     OTHER_best_Ag_j_class = OTHER_Ag_j_class_aff(1:OTHER_train_memory_ratio); 
% %                     OTHER_avg_aff = mean([OTHER_best_Ag_j_class]);

% %                     if(OTHER_avg_aff > avg_aff)
%                     if(OTHER_Ag_j_class_aff(1) > avg_aff)
                    max_other_Ag = max(OTHER_Ag_j_class_f_j);
                    
                    if(max_other_Ag > min_aff)
% % %                         total_good_neighboring_Ag = length(find ( [best_Ag_j_class] > max_other_Ag));
% % %                         if ( total_good_neighboring_Ag < min_neighborhood_size) 
% %                             temp = find(OTHER_Ag_j_class_f_j > (avg_aff - stdev) &&  OTHER_Ag_j_class_f_j < (avg_aff - 2* stdev));
%                           

                            bexist = true; % temporary breaker.... % Bad memory cell
                            break;
%                         else
%                             avg_aff = mean([best_Ag_j_class(1:total_good_neighboring_Ag)]);
%                             bexist = false;
%                             break;
%                         end
                    end
                end
            end
            
            

            if(~bexist && avg_aff > min_class_aff)
%                 mem_class_size = size(Ag_j_class_mapped_NRange,2);
%                 mem_set = reshape([Ab(Ag_j_class_mapped_NRange).receptor], [L mem_class_size])';
%                 str2 = num2str(Ab_mutated_best);       
% 
%                 bexist = false;
%                 for(i = 1:mem_class_size)
%                     str1 = num2str(mem_set(i,:));
%                     if(strcmp(str1(1,:) , str2(1,:)))
%                         bexist = true;
%                     end
%                 end



%                 if(~bexist)
                    Ab(min_class_loc_in_NRange).receptor = Ab_mutated_best;
                    Ab(min_class_loc_in_NRange).affinity = avg_aff;
                    f_j(min_class_loc_in_NRange) = avg_aff; % NOT NEEDED
        %             Ab(min_class_loc_in_NRange).affinity = best_aff; % NOT NEEDED
        %             since you are using f() for current afffinity

                % remove these antigens                 
                %<<  
                Ag(Ag_j_class_locations(Ag_j_class_position(1:size_best_Ag_j_class))) = [];
                j = j - length(find([Ag_j_class_locations(Ag_j_class_position(1:size_best_Ag_j_class))] <= j));                
                M = M - size_best_Ag_j_class;      

                if(M <= 0 || j>=M)
                   break; 
                end

                [class_info] = update_class_data(Ag, init_class_info);                
                %>>
                
                break;
            end

            % Use CLONALCLAS technique
            for ( i = 1:r)
                Ab(i).receptor = Ab_mutated(best_positions(i),:);
            end            
            
        end
        

        
        j = j + 1;
        
        Abd = generate_random_numbers(d, L, min_val, max_val, data_type);        
        
        % step 8
        [aff position] = sort(f_j(1:r),'ascend');
        for re = 1:d % these are non memory cells
            Ab(position(re)).receptor = Abd(re,:);
        end
        
        clear position;       
    end
    clear Abn
    
    tend = toc;
    fprintf('Gen %d elapsed time is %1.2f seconds\n',t, tend);

% % % % % % % %     fprintf('class: %d\n',init_class_info(1).number);
% % % % % % % %     [sort([Ab([init_class_info(1).Antibody_NRange]).affinity])]'
% % % % % % % %     
% % % % % % % %     fprintf('class: %d\n',init_class_info(2).number);
% % % % % % % %     [sort([Ab([init_class_info(2).Antibody_NRange]).affinity])]'
end

% output
% matured antibodies

%display according to your requirement
% % %<<
% %     for(i = 1:m)
% %         matured_antibodies{i} = reshape([Ab(r+i).receptor],[10 12])';
% %     end
% % %>>

% % % mem = [];
% % % for(i = 1:m)
% % %     mem(i,:) = [Ab(r+i).receptor Ab(r+i).class]; % --> not needed --> Ab(r+i).affinity/max_aff_XX*100];
% % % end

mem = [];
count = 1;
for(i = 1:m)
    if(Ab(r+i).affinity > min_memory_affinity)
        mem(count,:) = [Ab(r+i).receptor Ab(r+i).class];
        count = count + 1;
    end
end

file_used = ['data\' file_used '_memory.out'];
file_out = fopen(file_used,'w');

% fprintf(file_out,file_format, reshape([Ab([r+1:r+m]).receptor],[m L+1]));
fprintf(file_out, '%d \n\n', size(mem,2));
fprintf(file_out,file_format, mem');

fclose(file_out);

out = [[Ab(:).affinity]' [Ab(:).class]'];

fprintf('Preidicted Class ---- Given Class\n');
fprintf('      %1.4f                    %d\n',out');

load gong;
% wavplay(y,Fs);

