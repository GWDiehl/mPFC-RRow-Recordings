function keys = evalKeysFile(keysFunc)

% Function that produces a keys file from the provided file name as a
% string. Keys file should be built as a matlab function such that it can
% just be run using eval, however eval scares me so I am putting it in a
% function to prevent any funny buisness.

if ~exist(keysFunc,'file')
    error('The file does not exist is the current directory. Something is wrong, Stop Now!')
end
keys = eval(keysFunc);