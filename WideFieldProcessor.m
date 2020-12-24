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
    
    % Reading methods & pre-processing
    methods (Access = public)
        
        function readStack(obj, filename)
            
            disp('Reading tif stack...')
            tic
            if nargin == 1
                obj.Stack = TIFFStack(obj.filename);
            else
                obj.Stack = TIFFStack(filename);
            end
            time = toc;
            disp([' Time reading is ', time,' seconds.']);
            
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
            obj.f0 = obj.avg_projection;
            
            obj.dff = 100 * (image_raw - obj.f0) ./ obj.f0 ;
            
        end
        
    end
    
    % Post-Processing methods
    methods (Access = public)
        
        function registerSession(obj, session_reference, mouse_reference)
            
            
            
        end
        
        function  maskSession(obj, mask_struct)
            
        end
        
    end
    
    % Graphic methods
    
    methods
        
        function showPhotoBleaching(obj)
            
            num_frames = size(obj.Stack, 3);
            
            figure;
            plot(obj.frame_F,'linewidth',2)
            F_fit = polyfit(1:num_frames, obj.frame_F, 1);
            PhotoBl = round((F_fit(1) * num_frames) / F_fit(2) * 100);
            ylim([0 max(obj.frame_F) * 1.25])
            xlabel('Frame #')
            ylabel('Raw fluorescence')
            set(gcf,'color',[1 1 1])
            set(gca, 'FontSize', 14)
            title(['Fluorescence timecourse, Photobleaching = ' ...
                num2str(PhotoBl) '%'], 'FontSize', 18)
            
        end
        
        function showDFFMovie(obj, varargin)
            
            parser = inputParser;
            addParameter(parser, 'DownSamp', 1)
            parse(parser,varargin{:})
            
            down_samp = parser.Results.DownSamp;
            num_frames = size(obj.Stack, 3);
            max_value = 0.35 * max(obj.dff,[],'all');
            min_value = min(obj.dff, [], 'all');
           
            fig = figure;
            
            set(fig, 'visible', 'off');
            set(fig,'Color',[0 0 0])
            set(fig,'Position',[100 100 525 525]);
            subplot('Position',[0.05 0.05 0.9 0.9])
            colormap('jet');
            set(fig, 'visible', 'on');
            
            for i = 1:num_frames
                
                if ~ishghandle(fig)
                    break
                end
                
                if rem(i,down_samp) == 0
                    curr_frame = obj.dff(:,:,i);
                    imagesc(curr_frame);
                    caxis([min_value max_value])
                    axis square
                    title(['Frame ' num2str(i) '/' num2str(num_frames)],...
                        'color',[1 1 1]);
                    
                    drawnow()
                end

            end
            
        end
        
        
    end
    
end