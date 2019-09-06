function [ data data_type data_separator dist_type mutation_type file_used file_format] = preprocess_data( bprocess, fold )
%UNTITLED1 Summary of this function goes here
%  Detailed explanation goes here




%% data preprocessing
%<<
    %step 1: describe data
    data_type = '%g';
    data_separator = ',';
    dist_type = 'euclidean';
    mutation_type = 'euclidean';
    
       
%    data_type = '%d';
%    data_separator = ' ';
%    dist_type = 'hamming';
%    mutation_type = 'bit flip';
    
    
    %step 2: get raw data
    if(bprocess)
        file_type = '*.data';
    else
        file_type = 'train_*.data'
    end
    [data file_used] = get_data(data_type,data_separator,'Raw Data File',file_type);
    data_size = size(data,1);
    total_features = size(data,2) - 1; 
    
    %file format code
    %<<.................................................
        file_format = '';
        for(i = 1:total_features)
            file_format = [file_format data_type data_separator];
        end
        file_format = [file_format data_type];
        file_format = [file_format '\n'];
    %>>................................................
    
    if(~bprocess)
        return
    end

    
    
    
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
    sorted_data = sortrows(data, total_features + 1); % sort according to class
    class_list = unique(data(:,total_features + 1));
    
    count = 0;
    for (i = [class_list]')
        count = count + 1;
        class_size(count) = length(find (data(:,total_features + 1) == i));        
    end       
    
    % Code for PERCENTAGE SPLIT
    %<<    
% % %     cur_pointer = 0;
% % %     train_data = [];
% % %     test_data = []; 
% % %     training_ratio = 0.8;
% % %     for (i = 1:count)
% % %         rand_order = cur_pointer + randperm(class_size(i));
% % %         cur_pointer = cur_pointer + class_size(i);
% % %         
% % %         training_size = round(class_size(i) * training_ratio);
% % %         testing_size = class_size(i) - training_size;
% % %         train_data = [train_data' sorted_data(rand_order(1:training_size),:)']';
% % %         test_data = [test_data' sorted_data(rand_order(training_size+1:class_size(i)),:)']';    
% % %     end
    %>>
    
    % Code for cross-validation many folds
    %<<.................................................................    
    test = cell(fold,1);
    train = cell(fold,1);    
    cur_pointer = 1;
    
    for ( j = 1:count)
        a = sorted_data(cur_pointer:cur_pointer+class_size(j)-1,:);
        cur_pointer = cur_pointer + class_size(j);        
        sz = size(a,1);
        ratio = round(sz/fold);
        from = 0;
        
        for (i = 1:fold)
            to = from+1;
            from = min(i*ratio,sz);
            test{i} = [test{i}' a(to:from,:)']'; 
            temp = [a(1:to-1,:)' a(from+1:sz,:)']';
            train{i} = [train{i}' temp']';
        end
    end    
    %>>................................................................
    




    %<<    ????
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
    % Continuation Code for PERCENTAGE SPLIT
    %<<
% % %         train_file = ['data\train_' file_used];
% % %         test_file = ['data\test_' file_used];
% % % 
% % %         file_format = '';
% % %         for(i = 1:total_features)
% % %             file_format = [file_format data_type data_separator];
% % %         end
% % %         file_format = [file_format data_type];
% % %         file_format = [file_format '\n'];
% % % 
% % %         train_out = fopen(train_file,'w');    
% % %         fprintf(train_out, '%d \n\n', total_features + 1);
% % %         fprintf(train_out,file_format, train_data');
% % %         fclose(train_out);
% % % 
% % %         test_out = fopen(test_file,'w');    
% % %         fprintf(test_out, '%d \n\n', total_features + 1);
% % %         fprintf(test_out,file_format, test_data');
% % %         fclose(test_out);
    %>>
    
    % Continuation Code for cross-validation many folds
    %<<.................................................................  
    for i = 1:fold  
        train_file = ['data\train_' num2str(i) 'of' num2str(fold) '_' file_used];
        test_file = ['data\test_' num2str(i) 'of' num2str(fold) '_' file_used];

        train_out = fopen(train_file,'w');    
        fprintf(train_out, '%d \n\n', total_features + 1);
        fprintf(train_out,file_format, train{i}');
        fclose(train_out);

        test_out = fopen(test_file,'w');    
        fprintf(test_out, '%d \n\n', total_features + 1);
        fprintf(test_out,file_format, test{i}');
        fclose(test_out);
    end
    %>>.................................................................
    
    if(bprocess)
        error('data saved successfully. Rerun to process training data');
    end
    
%>>