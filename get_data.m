function [ data, filename] = get_data( data_type,  separator , msg, file_type)
%UNTITLED1 Summary of this function goes here
%  Detailed explanation goes here

[filename path] = uigetfile(file_type, msg);

if(path == 0)
    error(['file not selected' ' :-(']);
    exit
end

filepath = [path filename]

[fid,msg] = fopen(filepath,'r');

if(fid<0)
    error(['Cannot open ' filename]); 
    exit
end

fprintf('file (%s) opened successfully\n', filename);

L = fscanf(fid, '%d', [1 1]); % Total features/dimension with class number


format = '';
for(i = 1:L-1)
    format = [format data_type separator];
end
format = [format data_type];


%% Write code to get the data of your type
%<<
    file_data = fscanf(fid, format, [L inf]);
    file_data = file_data';
%>>

fclose(fid);
Antigen = file_data;

data = Antigen;
