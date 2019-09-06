function [vals, f] = get_frequencies(array)


if (size(array,1) ~= 1)
    array = array';
    if(size(array,1) ~= 1)
        error('Cannot get frequency for more then one dimension\n');
    end
end
array = sort(array);
[vals loc] = unique(array);
loc = [0 loc];
f = diff(loc);

end