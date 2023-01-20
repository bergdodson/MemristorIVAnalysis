function excelfilemaker(path, Data)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%Designating the eXcel filename, excel sheet name, save directory
pathSplit = strfind(path,'\');
saveSheetName = path(pathSplit(end)+1:end);
saveDirectory = path(1:pathSplit(end)); %strcat(filepattern(1:end-5),'D3 2-3_Data.xlsx');
saveFilename = strcat(path(pathSplit(end-1)+1:pathSplit(end)-1), '.xlsx');

measurementdata = struct;
dataFields = fields(Data);
%Prealocating memory for the data excel table
for x = 1:numel(dataFields)
    measurementdata.(dataFields{x}) = string(zeros(numel(Data),1));
end

%filling the preallocated memory variables created above 
%Some of the entries in Data are null values, some have multiple entries.
%to accomodate for each of the possibilities (no entry, 1 entry, multiple
%entries) the if statements found below were used.
for x = 1:numel(Data)
    %Filling in the data fields for each measurement
    for y = 1:numel(dataFields)
        measurementdata.(dataFields{y})(x) = exceldatatype(Data(x), dataFields{y});
    end    
end

%determining what the fields are that I want to save
fields2save = ["Filename","Vstart","Vend","SweepType","VsetRawIndex",...
    "VsetRaw","VsetDiff","VsetStatus","Vset","VresetRawIndex",...
    "VresetRaw","VresetDiff","VresetStatus","Vreset", "R100F","R100R",...
    "R100ratio","Vmax","Vmin","VsetInd","VresetInd"];

%Seeing which of the desired fields from above the measurement file actually contains
selectedData = find(arrayfun(@(j) ~isempty(find(string(fields2save) == string(j), 1)), dataFields));
variableHeader = dataFields(selectedData);

%turning all the data into a matrix
%preallocating memory
excelarray = string(zeros(numel(Data),numel(selectedData)));
for x = 1:numel(selectedData)
    %passing data to the matrix
        excelarray(:,x) = measurementdata.(dataFields{selectedData(x)});
end

%Turning the data into a table for excel
sampledata = array2table(excelarray, 'VariableNames', variableHeader);

%writing the table to an excel sheet
writetable(sampledata,strcat(saveDirectory,saveFilename),'Sheet',saveSheetName,'Range','A1'); %savefilename,'Sheet','Test','Range','A1')
  
end

