classdef WideFieldProcessor < handle
    % WidefieldProcessor Reads fluorescence activity from multi-tif files.
    %
    % Developed by SAM on 12/23/2020. Based on code developed by MJG and KKS.
    % Last updated 1/11/2020.
    
    properties (Access = public)
        Options         struct
        Stack
        filename(1,:)   char
        dff
        avg_projection
        f0
        frame_F
        mask            struct
        tform
        recovery_vec
        moving_time     logical
    end
    
    properties (Access = public)
        isRegistered    logical
        isMasked        logical
        isZipped        logical
    end
    
    % Initialization methods
    methods (Access = public)
        
        function obj = WideFieldProcessor(filename, varargin)
            
            parser = inputParser;
            addParameter(parser, 'Mode','auto')
            addParameter(parser, 'Frames', 'off')
            parse(parser ,varargin{:})
            
            obj.isRegistered = false;
            obj.isMasked = false;
            obj.isZipped = false;
            obj.tform = [];
            obj.moving_time = [];
        
            if nargin == 0
                obj.Options.Mode = 'manual';
                return
            else
                obj.getFile(filename);
            end
            
            transferOptions(obj, parser)
            
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
            
            if nargin == 1
                obj.Stack = TIFFStack(obj.filename);
            else
                obj.Stack = TIFFStack(filename);
            end
            
            
        end
        
        function computeProperties(obj)

            image_raw = single(obj.Stack(:, :, :));
            obj.avg_projection = mean(image_raw, 3);
            obj.frame_F = squeeze(mean(image_raw, [1 2]));
            
            if ~isempty(obj.moving_time)
                image_raw(:, :, obj.moving_time) = [];
                obj.f0 = mean(image_raw, 3);
            else
                obj.f0 = median(image_raw, 3);
            end

        end
        
        function computeDFF(obj)

            if ~isempty(obj.moving_time)
                image_raw = single(obj.Stack(:, :, obj.moving_time));
            else
                image_raw = single(obj.Stack(:, :, :));
            end
            
            if obj.isRegistered
                image_raw = imwarp(image_raw, obj.tform, 'OutputView', ...
                imref2d(size(image_raw)));
            end
            
            obj.dff = 100 * (image_raw - obj.f0) ./ obj.f0;
            
            if obj.isMasked
                obj.isMasked = false;
                obj.maskSession
            end
                
        end
            
    end
    
    % Post-Processing methods
    methods (Access = public)
        
        function registerDFF(obj)
            
            if obj.isRegistered == true
                disp('DFF has already been registered')
                return
            end
            
            if isempty(obj.tform)
                [moving_reference, fixed_reference] = getReferences(obj);
                obj.tform = obj.registerSession(moving_reference, fixed_reference);
            end
            
            obj.avg_projection = imwarp(obj.avg_projection, obj.tform, 'OutputView', ...
                imref2d(size(obj.avg_projection)));
            obj.f0 = imwarp(obj.f0, obj.tform, 'OutputView', ...
                imref2d(size(obj.f0)));
           
            obj.isRegistered = true;
            
            try
                obj.dff = imwarp(obj.dff, obj.tform, 'OutputView', ...
                    imref2d(size(obj.dff)));
            catch
                obj.isRegistered = false;
            end
            
        end
        
        function getMask(obj)
            % Gets the mask automatically
            
            session = strsplit(obj.filename, {'.', '_', '-'});
            mouse = session{1};
            mask_file = strcat(mouse,'_mask.mat');
            
            isSuccess = false;
            while ~isSuccess
                try
                    obj.mask = importdata(mask_file);
                    isSuccess = true;
                catch
                    [mask_file, ~] = uigetfile('*.mat');
                    
                    if isequal(obj.filename,0)
                        error('User selected Cancel')
                    end
                end
            end
            
        end
        
        function maskSession(obj)
            % It only masks the hemispheres
            
            if obj.isMasked
                disp('Session has already been masked')
                return
            end
            
            if ~obj.isRegistered
                disp('Session has to be registered first')
                return
            end
            
            if isempty(obj.mask)
                obj.getMask
            end
             
            mask_mat = obj.mask.ALL;
            mask_mat(~obj.mask.Window) = 0;
            
            if obj.isZipped
                obj.unzipDFF
            end
            
            mask_mat_DFF = repmat(mask_mat, [1 1 size(obj.dff,3)]);
            
            obj.f0(~mask_mat) = 0;
            obj.avg_projection(~mask_mat) = 0;
            
            try
                obj.dff(~mask_mat_DFF) = 0;
                obj.isMasked = true;
            catch
                obj.isMasked = false;
            end
            
        end
        
        function zipDFF(obj)
            
            if obj.isZipped
                disp('DFF already zipped')
                return
            end
            
            if isempty(obj.dff)
                disp('DFF not available')
                return
            end
            
            [obj.dff, obj.recovery_vec] = obj.zipMat(obj.dff);
            obj.isZipped = true;
           
        end
        
        function unzipDFF(obj)
            
            if ~obj.isZipped
                disp('DFF already not-zipped')
                return
            end
            
            if isempty(obj.dff)
                disp('DFF not available')
                return
            end
            
            [obj.dff, ~] = obj.zipMat(obj.dff, obj.recovery_vec);
            obj.isZipped = false;
            
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
   
    % Methods for general use
    methods (Access = public, Static)
        
        function tform = registerSession(moving_reference, fixed_reference)
            % registerSession Registration of the session reference with
            % the overall mouse reference.
            %
            % Input must be of type uint16. If it is double, it will be
            % enough to call it as 'uint16(moving_image)'
            % Output is the transformation tensor. In order to transform
            % any matrix or image, it can be used as:
            % imwarp(moving_image, tform, 'OutputView', imref2d(size(moving_image)))
           
            [moving_pts, fixed_pts] = cpselect(moving_reference, ...
                fixed_reference, 'Wait',true);
            
            if ischar(moving_reference)
                moving_reference = imread(moving_reference);
                fixed_reference = imread(fixed_reference);
            end
            
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
            
            fixed_image = imwarp(moving_reference, tform, 'OutputView',...
                imref2d(size(moving_reference) ));
            
            try
                [fig_opt, ax_opt, t_opt] = figure_properties();
                
                figure(fig_opt{:}, 'Position', [200 200 1200 600])
            catch
                figure('Position', [200 200 1200 600]');
            end
            
            subplot(1, 2, 1)
            imshowpair(moving_reference, fixed_reference)
            set(gca, ax_opt{:}),
            title('Original Reference', t_opt{:});
            
            subplot(1, 2, 2)
            imshowpair(fixed_image, fixed_reference)
            set(gca, ax_opt{:}),
            title('Processed', t_opt{:});
            
        end
        
        function [mat_new, recovery_vec] = zipMat(mat_old, recovery_vec)
            % The goal of this subroutine is to reduce the size of large matrices that
            % have been masked with zeros. It removes the zeros so that it can be saved
            % (along with the recovery vector). If the recovery vector is introduced,
            % then it takes the matrix back to its normal shape. The function assumes
            % the input matrices to be either [x_pixels y_pixels n_frames] or [n_pixels
            % n_frames].
            %
            % Developed by Santi on 29/12/2020
            % Last updated 3/1/2020
            
            if nargin == 1
                
                n_pixels = size(mat_old, 1) * size(mat_old,2);
                n_frames = size(mat_old, 3);
                
                mat_old = reshape(mat_old, [n_pixels, n_frames]);
                recovery_vec = mat_old(:,1) ~= 0;
                
                mat_new = mat_old;
                mat_new(~recovery_vec, :) = [];
                
            elseif nargin == 2
                
                n_pixels = length(recovery_vec);
                n_frames = size(mat_old, 2);
                
                mat_new = zeros(n_pixels, n_frames, 'single');
                
                k = 0;
                for i = 1:n_pixels
                    
                    if recovery_vec(i)
                        k = k + 1;
                        mat_new(i,:) = mat_old(k,:);
                    end
                    
                end
                
                mat_new = reshape(mat_new, [sqrt(n_pixels) sqrt(n_pixels) ...
                    n_frames]);
                
            end
            
        end
        
    end
    
    % Graphic methods
    methods (Access = public)
        
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
        
        function showDFFMovie(obj, DFF, varargin)
            
            if nargin == 1
                
                DFF = obj.dff;
                
                if obj.isZipped
                    DFF = obj.zipMat(DFF, obj.recovery_vec);
                end
                
            end
            
            min_val = prctile(obj.dff(:), 1);
            max_val = prctile(obj.dff(:), 99);
            
            parser = inputParser;
            addParameter(parser, 'DownSamp', 1)
            addParameter(parser, 'Mask', 'window')
            parse(parser,varargin{:})
            
            down_samp = parser.Results.DownSamp;
            num_frames = size(DFF, 3);
            
            if obj.isRegistered
                switch parser.Results.Mask
                    
                    case 'off'
                        mask_region = ones(400);
                        
                    case 'window'
                        obj.getMask
                        mask_region = obj.mask.ALL;
                        mask_region(~obj.mask.('Window')) = 0;
                        
                    otherwise
                        obj.getMask
                        window_mask = obj.mask.('Window');
                        mask_region = obj.mask.(parser.Results.Mask);
                        mask_region(~window_mask) = 0;
                end
                
            else
                disp('Mask not displayed because the session is not registered')
                mask_region = ones(400);
            end
                            
            fig = figure;
            
            set(fig, 'visible', 'off');
            set(fig,'Color',[0 0 0])
            set(fig,'Position',[100 100 525 525]);
            
            subplot('Position',[0.05 0.05 0.9 0.9])
            cmap = parula(500);
            colormap(cmap)
            set(fig, 'visible', 'on');
            
            for i = 1:num_frames
                
                if ~ishghandle(fig)
                    break
                end
                
                if rem(i,down_samp) == 0
                    curr_frame = DFF(:,:,i);
                    
                    imagesc(curr_frame, 'AlphaData', mask_region);
                    
                    caxis([min_val max_val])
                    axis square
                    title(['Frame ' num2str(i) '/' num2str(num_frames)],...
                        'color',[1 1 1]);
                    set(gca,'color',0*[1 1 1]);
                    
                    
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