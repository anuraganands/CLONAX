function [ C_mutated ] = hypermutate( aff, ro, C, L, type, varargin)
    % TODO........
    % mutation can occur in different manner depending on the input
    % possible mutation types are
    % 'bit flip'  -   when a binary bit is inverted
    % 'val swap'  -   when two randomly picked values are swapped  ... NO
    % LONGER USED
    % 'euclidean' - variable multiplied by [0 1] range of values.

    norm_data = norm_scale01(aff); % normalize data over the interval [0,1]
    mutation_rate = exp(-ro*norm_data); % rate between (0,1]
    n = size(C,1); % clones of n unique antibodies
    
    if(strcmp(type,'euclidean') == 1)
        if(nargin ~= 7)
            error('min and max value not specified');
        end
        min_val = varargin{1};
        max_val = varargin{2};
    else
        if(nargin == 5)
            min_val = []; % Value not needed
            max_val = []; % Value not needed
        end
    end
    
    for(cl = 1:n) % clones of n unique antibodies
    %     clone_size = size(C{j,1}{cl,1},1); %just getting row size which is clone size
        clone_size = size(C{cl,1},1);
        alpha = mutation_rate(cl);
    %             data = C{j,1}{cl,1};

        % mutation rate is designed in such a way so that it is
        % inversely proportional to affinity. i.e. Higher affinity has
        % low mutation rate and vice versa.
        mutation_freq = ceil(alpha*L); %alpha is mutation rate and L is lenght of Antibody

        % mutate alpha times to given antibody
        % get all antibodies
        if(mutation_freq > 0)        
            for( a = 1:clone_size)
                for(i = 1:mutation_freq)
                    swap1 = round( random_locator(1, L) );

                    if(strcmp(type,'bit flip') == 1)                        
                        bit = C{cl,1}(a,swap1);
                        if(bit == 0)
                            bit = 1;
                        else
                            bit = 0;
                        end                        
                        C{cl,1}(a,swap1) = bit; %flip the bit                        
                    elseif(strcmp(type,'val swap') == 1)
                        swap2 = round( random_locator(1,L) );
                        temp = C{cl,1}(a,swap1);
                        C{cl,1}(a,swap1) = C{cl,1}(a,swap2);
                        C{cl,1}(a,swap2) = temp;
                    elseif(strcmp(type,'euclidean') == 1)
                        eu_val = C{cl,1}(a,swap1);
                        C{cl,1}(a,swap1) = min_val + mod(rand(1)*1.0*(max_val - min_val) + eu_val, max_val - min_val); % only 10% i.e. 0.1 of the range is added at a time           
%                         C{cl,1}(a,swap1) = min_val + rand(1)*(max_val - min_val); % only 10% i.e. 0.1 of the range is added at a time           
                    else
                        error('Incorrect argument value');
                    end
                end
            end
        end
    end

    C_mutated = C;
    
end


function r = random_locator(min_val, max_val)
    r = min_val +(max_val - min_val)* rand(1);
end

