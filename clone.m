function [ C, array_clone_size ] = clone( Abnj, beta, N )
% Clone given n antibodies

n = size(Abnj,1);

C = cell(n,1);

for(i = 1: n)
    clone_size = round(beta*N/i); % different clone sizes for each of antibodies depending on affinity strength.
    array_clone_size(i) = clone_size;
    for(j = 1: clone_size)
       C{i,1}(j,:) = Abnj(i,:); 
    end
end

