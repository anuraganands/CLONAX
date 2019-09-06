
function [ struct_out ] = clonal_struct( struct_desc )
%UNTITLED1 Summary of this function goes here
%  Detailed explanation goes here

if ischar(struct_desc), 
    switch struct_desc
        case 'Antigen',
            struct_out = struct(    'type','Antigen', ...
                                    'epitope',[], ... %epitope
                                    'class',-1);
        case 'Antibody',
            struct_out = struct(    'type','Anitibody', ...
                                    'receptor',[], ... %receptor
                                    'class',0, ...
                                    'affinity',-1, ...
                                    'isMemoryCell',false);
        case 'class_info',
            struct_out = struct(    'type','class_info', ...
                                    'number',-1, ...
                                    'size',-1, ...
                                    'Antibody_NRange', [0 0],...
                                    'locations',[]);
        otherwise
            error('Unrecognized argument for struct_desc');
            struct_out = [];
            return;
    end  
else
    error('Incorrect type for struct_desc');
    struct_out = [];
    return;
end