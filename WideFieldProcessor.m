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
    
    properties (Access = public)
        isRegistered logical
        isMasked logical
    end
    
    % Initialization methods
    methods (Access = public)
        
        function obj = WideFieldProcessor(filename, varargin)
            
            parser = inputParser;
            addParameter(parser, 'Mode','auto')
            addParameter(parser, 'Frames', 'off')
            parse(parser ,varargin{:})
            transferOptions(obj, parser)
            
            obj.isRegistered = false;
            
            
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
    
    %% Post-Processing methods
    methods (Access = public)
        
        function registerDFF(obj)
            
            if obj.isRegistered == true
                disp('DFF has already been registered')
                return
            end
            
            [moving_reference, fixed_reference] = getReferences(obj);
            tform = registerSession(moving_reference, fixed_reference);
            
            obj.dff = imwarp(obj.dff, tform, 'OutputView', ...
                imref2d(size(obj.dff)));
            
            obj.isRegistered = true;
            
        end
        
        function maskDFF(obj)
            
            if obj.isMasked == true
                disp('DFF has already been masked')
                return
            end
            
        end
        
    end
        
       
    methods (Access = private)
        
        function [moving_reference, fixed_reference] = getReferences(obj)
            % findReferences Looks for both references. If the absolute
            % reference (based on the filename) is not found, it prompts
            % the user to specify one. If the session reference is not
            % found, it uses the average projection as reference. 
           
            session = strsplit(obj.filename, {'.', '_', '-'});
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
                            moving_reference = uint16(obj.avg_projection);
                            break
                            
                        otherwise
                            warning('Only 1 or 2 options are accepted');
 
                    end

                end
            end

        end
        
        
    end
    
    
    methods (Access = public, Static)
        
        function tform = registerSession(moving_image, fixed_reference)
            % registerSession Registration of the session reference with
            % the overall mouse reference.
            %
            % Input must be of type uint16. If it is double, it will be
            % enough to call it as 'uint16(moving_image)'
            % Output is the transformation tensor. In order to transform
            % any matrix or image, it can be used as:
            % imwarp(moving_image, tform, 'OutputView', imref2d(size(moving_image)))
            
            [moving_pts, fixed_pts] = cpselect(moving_reference, fixed_reference,...
                'Wait',true);
            
            tform = fitgeotrans(moving_pts, fixed_pts, ...
                'nonreflectivesimilarity');
            
            u = [0 1];
            v = [0 0];
            [x, y] = transformPointsForward(tform, u, v);
            dx = x(2) - x(1);
            dy = y(2) - y(1);
            
            angle = (180/pi) * atan2(dy, dx);
            scale = 1 / sqrt(dx^2 + dy^2);
            
            disp([' Rotation angle is ', num2str(angle), ' degrees.'])
            disp([' Scale is ', num2str(scale), '.'])
            
            fixed_image = imwarp(moving_image, tform, 'OutputView',...
                imref2d(size(moving_image) ));
            
            try
                [fig_opt, ax_opt, t_opt] = figure_properties();
                
                figure(fig_opt{:}, 'Position', [200 200 1200 600])
            catch
                figure('Position', [200 200 1200 600]');
            end
            
            subplot(1, 2, 1)
            imshowpair(moving_image, fixed_reference)
            set(gca, ax_opt{:}),
            title('Original Reference', t_opt{:});
            
            subplot(1, 2, 2)
            imshowpair(fixed_image, fixed_reference)
            set(gca, ax_opt{:}),
            title('Processed', t_opt{:});
            
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
    
    % Other methods
    
    methods (Access = private, Static)
        
        function isCancelled(file)
            
            if isequal(file, 0)
                error('User selected Cancel')
            end
            
        end
        
    end
    
end