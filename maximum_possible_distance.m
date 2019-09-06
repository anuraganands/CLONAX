function [ max_dist ] = maximum_possible_distance( data )
% get minimum and maximum for each feature to derive lowest and highest
% possible value.


M = size(data,1); %rows
L = size(data,2); % cols => features/dimension
min_array = [];
max_array = [];
for(i = 1:L)
   min_array(i) = min(data(:,i)); 
   max_array(i) = max(data(:,i));
end

max_dist = norm(max_array - min_array);






% % % %ONLY IF DATA SET IS NORMALIZED WITH [0 1] ********************************
% % % sq = 0;
% % % for(i = 1:L)
% % %     sq = sq + 1*1;
% % % end
% % % max_dist = sqrt(sq);
