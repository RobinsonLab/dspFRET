function [] = writeStruct(filePtr,s)
% INPUT
% filePtr, file pointer. 1 for stdio 
% s, structure consisting of name, value
% OUTPUT
% none
fields = fieldnames(s);
for i=1:length(fields)
    fprintf(filePtr, '%s %d \n',fields{i}, getfield(s,fields{i}));
end
end