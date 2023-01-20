function output = exceldatatype(Data,datafield)
%Function: exceldatatype.m
%Description: Determines how many elements are in a given measurement
%datafield (0, 1, or many elements). The function then arranges the data
%appropriately to be writting into an excel table

if ~isempty(Data.(datafield))
    if numel(Data.(datafield)) == 1
        output = Data.(datafield);
    else
        output = strjoin(string(Data.(datafield)),', ');
    end
else
    output = ""; %strjoin(string([]))
end
end

