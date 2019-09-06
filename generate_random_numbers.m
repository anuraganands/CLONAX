function [ X ] = generate_random_numbers( row, col, min_val, max_val, data_type )
%UNTITLED1 Summary of this function goes here
%  Detailed explanation goes here



if(strcmp(data_type,'%d') == 1)
    X = min_val +(max_val-min_val)*randint(row,col); % binary numbers
else
    X = min_val +(max_val-min_val)*rand(row,col);
end