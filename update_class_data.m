function [ uclass_info varargout] = update_class_data( Ag, class_info, Ab, m, r, M)
% Ag - current Antigens


class_list = unique([Ag(:).class]); % class numbers

bupdate_mem = false;
if(nargin > 2)
    bupdate_mem = true;
else
    m = -1;
    r = -1;
    M = -1;    
end

N = m + r;
counter = 1;
Nidx = r;

if(bupdate_mem)
    for(i = [class_list])
        class_info(counter).type = 'class_info';
        class_info(counter).number = i;
        class_info(counter).locations = find([Ag(:).class] == i);
        class_info(counter).size = size([class_info(counter).locations],2);
        
        memory_class_size = round(class_info(counter).size/M*m);
        class_info(counter).Antibody_NRange = [Nidx+1  :  Nidx+memory_class_size];        
        Nidx = Nidx + memory_class_size;

        counter = counter + 1;
    end  
    % Update last one so that rounding off doesn't cause the last memory block
    % to be unassigned.

    class_info(counter - 1).Antibody_NRange = [class_info(counter - 1).Antibody_NRange(1) : N];

    %assign class numbers to Ab
    counter = 1;
    for(i = [class_list]) 
        for(j = [class_info(counter).Antibody_NRange])
            Ab( j ).class = class_info(counter).number;
        end
        counter = counter + 1;
    end
    
    uclass_info = class_info;
    varargout{1} = Ab;
    varargout{2} = class_list;
else
    for(i = [class_list]) %it contains actual class number
        cf(counter).type = 'class_info';
        cf(counter).number = i;
        cf(counter).locations = find([Ag(:).class] == i);
        cf(counter).size = size([cf(counter).locations],2);
               
        cf(counter).Antibody_NRange = [class_info( find([class_info(:).number] == i)).Antibody_NRange];  %????????????? it will not change...     

        counter = counter + 1;
    end  
    % Update last one so that rounding off doesn't cause the last memory block
    % to be unassigned.  
    uclass_info = cf;
end

