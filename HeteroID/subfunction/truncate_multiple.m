function [y_trim, A_trim] = truncate_multiple(y, A, grouping, partition, cut_coefficient)
% Use cut_coefficient (in percentage) data of each condition to do the
% identification

%%
% clear all
% cut_coefficient = .4;

A_trim = A;
y_trim = y;
 

for i = 1:length(partition)
    length_grouping(i) = length(grouping{i});
end
[length_grouping_max_val,length_grouping_max_position] = max(length_grouping);


for j = 1:length_grouping_max_val
    condition_position(j) = grouping{length_grouping_max_position}(j);
    condition_position_nnz{j} = find(A(:,condition_position(j)));
    condition_length(j) = length(condition_position_nnz{j});
    condition_length_cut(j) = round(cut_coefficient*condition_length(j));
end

condition_length_cut_grouping = [];
condition_length_initial(1) = 1;
cumsum_condition_length = cumsum(condition_length) ;
clear j;

for j = 1:length_grouping_max_val
    condition_length_initial(j+1) =  cumsum_condition_length(j) + 1;
    condition_length_cut_grouping_j = ...
        condition_length_initial(j)+condition_length_cut(j): condition_length_initial(j+1)-1;
    condition_length_cut_grouping = [condition_length_cut_grouping, ...
        condition_length_cut_grouping_j ];
end


A_trim(condition_length_cut_grouping, :) = [];
y_trim(condition_length_cut_grouping, :) = [];