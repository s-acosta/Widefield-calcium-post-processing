% Test Widefield Processor Class


%% Initialazing test

% Check if it is able to read and get files
filename = 'WTR039_ref.tif';

WF_obj = WideFieldProcessor(filename, 'Frames', 'off');

