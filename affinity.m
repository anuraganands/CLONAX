function [af] = affinity(Ab, Ag_j, distance_type, varargin)
    % Ab is antibodies of size N = m + r. m => memory cells, r => common
    % antibody cells

    % Ag_j = One particular antigen

    L = size(Ab,2);
    N = size(Ab,1); %Antibody population size

    if(nargin > 5)
        error('Wrong number of input arguments') 
    end
    
    bnormalize = false;
    if(nargin >= 4)
        bnormalize = true;
        max_dist = varargin{1};
    end
    
    if(nargin == 5) % Currenty it can be used in testing phase of classification.
        ignore_error = varargin{2};
    else 
        ignore_error = 1; % Do not ignore anything...
    end
    
    if(strcmp(distance_type,'hamming')==1) % find Hamming distance (HD)
        for(n = 1: N) % for each antibody
            HD = 0;
            for e = 1: L % length of each antigen
                if(Ag_j(e) == Ab(n,e)) % check each epitope
                    HD = HD + 1;     
                else
                    ; % HD = HD + 0;
                end                
            end
            af(n) = HD; %higher hamming dist => lower affiniy
            if(bnormalize)
                af(n) = af(n)/max_dist; % normalized distance in range [0,1]
            end
        end
    elseif(strcmp(distance_type,'euclidean') == 1)
        for(n = 1:N) % for each antibody
            temp_Abn = Ab(n,:);
            temp_Ag_j = Ag_j;
            
% % %             if(ignore_error ~= 1)
% % %                 max_err = L - round(L*ignore_error);                
% % %                 [vals err_idx] = sort(abs(Ag_j-Ab(n,:)),'descend');
% % %                 % remove error features - i.e. slots with very high
% % %                 % difference
% % %                 temp_Ag_j(err_idx(1:max_err)) = [];
% % %                 temp_Abn(err_idx(1:max_err)) = [];
% % %             end 
            
            af(n) =  max_dist - norm(temp_Ag_j - temp_Abn); %higher euclidean dist => lower affiniy
            if(bnormalize)
                af(n) = af(n)/max_dist; % normalized distance in range [0,1]
            end
        end
    else
        error('Incorrect argument for distance_type');
    end
end