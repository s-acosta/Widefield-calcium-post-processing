% test_findReferences

average_projection = ones(3);
% filename = 'SAM000_201114'; 
filename = 'SAM000_000000';

session = strsplit(filename, {'.', '_', '-'});
mouse = session{1};
fixed_reference_file = strcat(mouse,'_ref.tif');

isSuccess = false;
while ~isSuccess
    try
        fixed_reference = imread(fixed_reference_file);
        isSuccess = true;
    catch ME
        warning(ME.message)
        disp('Please choose manually')
        [fixed_reference_file, ~] = uigetfile({'*.tif';'*.tiff'});
        WideFieldProcessor.isCancelled(fixed_reference_file)
    end
end

date = session{2};
moving_reference_file = strcat(mouse,'_',date,'_ref.tif');

isSuccess = false;
while ~isSuccess
    try
        moving_reference = imread(moving_reference_file);
        isSuccess = true;
    catch ME
        warning(ME.message)
        prompt = 'Choose file (1) or average projection (2)? ';
        answer = input(prompt);
        
        switch answer
            
            case 1
                [moving_reference_file, ~] = uigetfile(...
                    {'*.tif';'*.tiff'});
                WideFieldProcessor.isCancelled(moving_reference_file);
                
            case 2
                moving_reference = uint16(average_projection);
                break
                
            otherwise
                warning('Only 1 or 2 options are accepted');
                
        end
        
    end
end
