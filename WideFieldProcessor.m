classdef WideFieldProcessor < handle
    % WidefieldProcessor Reads fluorescence activity from multi-tif files.
    %
    % Developed by SAM on 12/23/2020. Based on code developed by MJG and KKS. 
    % Last updated 12/23/2020.
    
    properties (Access = public)
        Options  struct
        Stack
        filename(1,:) char
        dff
        avg_projection
        f0
        frame_F
    end
    
    % Initialization methods
    methods (Access = public)
        
        function obj = WideFieldProcessor(filename, varargin)
            
            parser = inputParser;
            addParameter(parser, 'Mode','auto')
            addParameter(parser, 'Frames', 'off')
            parse(parser ,varargin{:})
            transferOptions(obj, parser)
            
            if nargin == 0
                obj.Options.Mode = 'manual';
                return
            else
                obj.getFile(filename);
            end
            
        end
        
        function getFile(obj, filename)
           
            if nargin == 2
                obj.filename = filename;
            else
                [obj.filename, ~] = uigetfile({'*.tif';'*.tiff'});
                
                if isequal(obj.filename,0)
                    error('User selected Cancel')
                end
                
            end
            
        end
        
        function setOption(obj, optionName, optionValue)

            assert(ischar(optionName), 'First argument has to character type')
            
            error_msg = 'Invalid value for that option';
            
            switch optionName
                
                case 'Frames'
                    assert((ischar(optionValue) && contains(optionValue, 'off')) | ...
                        (isnumeric(optionValue) && ndims(optionValue) == 2), ...
                        error_msg);
                
                case 'Mode'
                    assert(ischar(optionValue), error_msg);
                    assert(contains(optionValue, 'auto') | ...
                        contains(optionValue, 'manual') , error_msg);
                    
                otherwise
                    error('That option does not exist')
                    
            end
            
            obj.Options.(optionName) = optionValue;
                
            
        end
            
    end
    
    methods (Access = private)
        
        function transferOptions(obj, parser)
            
            obj.Options = struct;
            obj.Options.Mode = parser.Results.Mode;
            obj.Options.Frames = parser.Results.Frames;
            
        end
        
    end
    
    % Reading methods
    methods (Access = public)
        
        function readStack(obj, filename)
                
                disp('Reading tif stack...')
                if nargin == 1
                    obj.Stack = TIFFStack(obj.filename);
                else
                    obj.Stack = TIFFStack(filename);
                end
                
        end
        
        function computeProperties(obj)
            
            disp('Computing basic properties...')
            
            switch obj.Options.Frames
                case 'off'
                    image_raw = single(obj.Stack(:, :, :));   
                otherwise
                    image_raw = single(obj.Stack(:, :, obj.Options.Mode));
            end
            
            obj.avg_projection = mean(image_raw, 3);
            obj.frame_F = squeeze(mean(image_raw, [1 2]));
            obj.f0 = median(image_raw, 3);
         
        end   
            
        
    end
    
    
        
        
        
    
    
end