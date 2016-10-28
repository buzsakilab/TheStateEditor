function StateEditor(baseName, inputData, supressGUI, makePortable)
%Just calls the properly named TheStateEditor function
% BW July 2015

if ~exist('baseName','var')
    baseName = [];
end

if ~exist('inputData','var')
    inputData = [];
end

if ~exist('supressGUI', 'var')
    supressGUI = 0;
end

if ~exist('makePortable', 'var')
    makePortable = 0;
end


TheStateEditor(baseName, inputData, supressGUI, makePortable)