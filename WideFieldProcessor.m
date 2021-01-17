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
        dff_deconvolved
        dff_sigma
        dff_baseline
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
        isThresholded   logical
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
            obj.isThresholded = false;
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
        
        function image_raw = readTiff(obj, frames)
            
           block_size = 500;
            
           if nargin == 1
               frames = true(1, size(obj.Stack,3));
           end
           
           n_frames = length(frames);
           
           if n_frames < block_size
               image_raw = single(obj.Stack(: ,:, frames));
           
           else
               
               image_raw = zeros(size(obj.Stack,1), size(obj.Stack,2), ...
                   n_frames, 'single');
               
               
               num_blocks = floor(n_frames / block_size);
               
               progressBar(0)
               for i = 1:num_blocks
                   
                   progressBar(i/num_blocks);
                   
                   startt = (i-1)*block_size + 1;
                   endd = i * block_size; 
                   
                   image_raw(:, :, startt : endd) = ...
                       single(obj.Stack(: ,:, frames(startt:endd)));
                   
               end
          
               if ~isequal(mod(n_frames, block_size) , 0)
                   image_raw(:, :, num_blocks*block_size + 1:end) = ...
                       single(obj.Stack(: ,:, ...
                       frames(num_blocks*block_size + 1:end)));
               end
                       
                       
           end
           
        end
        
        function computeProperties(obj)

            image_raw = obj.readTiff;
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

            image_raw = obj.readTiff;
            
            if obj.isRegistered
                image_raw = imwarp(image_raw, obj.tform, 'OutputView', ...
                imref2d(size(image_raw)));
            end
            
            obj.dff = 100 * (image_raw - obj.f0) ./ obj.f0;
            
            if ~isempty(obj.moving_time)
                obj.dff_sigma = std(obj.dff(:,:,~obj.moving_time), 1, 3);
                obj.dff_baseline = mean(obj.dff(:,:,~obj.moving_time), 3);
            else
                obj.dff_sigma = std(obj.dff, 3);
                obj.dff_baseline = mean(obj.dff, 3);
            end
            
        end
        
        function cropDFF(obj)
            
            if ~isempty(obj.moving_time)
                obj.dff = obj.dff(:, :, obj.moving_time);
            else
                disp('Cant be cropped because there is no moving_time')
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
        
        function thresholdDFF(obj)
            %thresholdDFF Thresholding of DFF traces
            % 
            % This method threshold the DFF traces by only considering the
            % ones that have maximum values higher than 2 times the
            % standard deviation of the baseline fluorescence. After
            % maximum values are computed, it keeps the transient until it
            % reaches 0.5 std deviation, then it sets to zero all the other
            % points until the next transient. Same procedure as Dombeck et
            % al, Neuron 56, 43-57, 2007.
            %
            % Developed by Santi on 1/15/2020
            % Last updated 1/15/2020
            
            if obj.isThresholded
                disp('It has already been thresholded')
                return
            end
            
            if isequal(ndims(obj.dff), 3)
                obj.zipDFF
            end
            
            valid_idxs = find(obj.recovery_vec);
            n_frames = size(obj.dff, 2);
            
            for j = 1:length(valid_idxs)
                
                std = obj.dff_sigma(valid_idxs(j));
                baseline = obj.dff_baseline(valid_idxs(j));
                isHigh = obj.dff(j,:) > (2*std + baseline);
                where_isHigh = find(isHigh);
                
                obj.dff(j, 1:where_isHigh(1)) = 0;
                
                for i = 1:length(where_isHigh)
                    
                    if i == length(where_isHigh)
                        temp = obj.dff(j, where_isHigh(i):n_frames);
                        obj.dff(j, where_isHigh(i):n_frames) = temp;
                        
                    elseif where_isHigh(i+1) - where_isHigh(i) == 1
                        continue
                        
                    else
                        temp = obj.dff(j, where_isHigh(i):where_isHigh(i+1));
                        where_ends = find(temp < 0.5 * (std + baseline), 1, 'first');
                        where_ends = where_ends + where_isHigh(i);
                        obj.dff(j, where_ends:where_isHigh(i+1)) = 0;
                        
                    end
                    
                end
                
            end
            
            obj.isThresholded = true;
          
           
        end
        
        function deconvolveDFF(obj)
            % Deconvolution by Lucy-Richardson method, as shown in Merav Stern, 
            % Eric Shea-Brown, Daniela Witten bioRxiv 2020.02.01.930040;
            % 
            % Developed by SAM on 1/15/2020
            % Last updated 1/15/2020
            
            if ~obj.isThresholded
                disp('First threshold the signal, works better')
                return
            end
            
            if ~obj.isZipped
                obj.zipDFF
            end
            
            obj.dff_deconvolved = zeros(size(obj.dff,1), size(obj.dff,2), ...
                'single');
            
            gamma = 0.95; % Decay parameter, taken for Clancy et al, 2019.
            p = 0:1:500; % Points
            conv_kernel = gamma.^p;
            kernel_for_lucy = [zeros(size(conv_kernel)) conv_kernel]';
            
            for i = 1:size(obj.dff,1)
                obj.dff_deconvolved(i,:) = deconvlucy(obj.dff(i,:)', ...
                    kernel_for_lucy);
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
                recovery_vec = mean(mat_old,2) ~= 0;
                
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
        
        function [traces, pixel_idx] = getTraces(mat, mask, n_pixels)
            % The goal of this subroutine is to extract traces from n
            % random pixels inside the mask matrix, which has to be a
            % logical with size pixels x pixels
            %
            % Developed by Santi on 1/14/2020
            % Last updated 1/14/2020
            
            mask_true = find(mask);
            
            pixel_idx = zeros(n_pixels, 2);
            traces = zeros(n_pixels, size(mat, 3));
            
            for i = 1:n_pixels
                
                element = randi(numel(mask_true));
                [pixel_idx(i,1), pixel_idx(i,2)] = ind2sub([size(mat,1), ...
                    size(mat,2)], mask_true(element));
                
               traces(i,:) = mat(pixel_idx(i,1), pixel_idx(i,2), :);
               
            end
            
            linear_idxs = sub2ind([size(mat,1) size(mat,2)], ...
                pixel_idx(:,1), pixel_idx(:,2));
            [~, order] = sort(linear_idxs, 'ascend');
            
            traces = traces(order, :);
            pixel_idx = pixel_idx(order, :);
            
        end
        
        function showRegionTraces(mat, mask, traces, pixel_idx)
            
            n_pixels = size(traces,1);
            pixel_mat = false(size(mat, 1), size(mat, 2));
            
            for i = 1:n_pixels
                pixel_mat(pixel_idx(i,1)-5 : pixel_idx(i,1) + 5, ...
                   pixel_idx(i,2) - 5 : pixel_idx(i,2) + 5) = true;
               pixel_mat(pixel_idx(i,1)-3 : pixel_idx(i,1) + 3, ...
                   pixel_idx(i,2) - 3 : pixel_idx(i,2) + 3) = false;
            end
                
            mask(pixel_mat == true) = false;
            
            figure('Color', [1 1 1], 'Position', [200 200 1200 700]);
            
            tiledlayout(n_pixels, 2*n_pixels);
            
            nexttile(1, [n_pixels n_pixels])
            imagesc(mean(mat,3), 'AlphaData', mask)
            ax_opt = {'XColor',[0 0 0], 'YColor',[0 0 0], ...
                'Colormap',colormap('parula'), 'DataAspectRatio', [1 1 1]};
            set(gca, ax_opt{:});
            set(gca,'Visible','off')
            
            ymax = max(traces(:)) + 5;
            ymin = min(traces(:)) + 5;
            
            for i = 1:n_pixels
                nexttile(n_pixels + 1 + n_pixels*2*(i-1), [1 n_pixels])
                plot(traces(i,:))
                ylim([ymin ymax]);
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
            
            min_val = 0;
            max_val = 30;
            
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
    
    methods (Access = public, Static)
        
        function showMap(mat, mask, lims)
            
            figure('Color', [1 1 1], 'Position', [200 200 600 600]);
            
            if nargin > 1
                imagesc(mat, 'AlphaData', mask)
            else
                imagesc(mat)
            end
            
            if nargin == 3
                caxis([lims(1) lims(2)])
            end
            
            ax_opt = {'XColor',[0 0 0], 'YColor',[0 0 0], ...
                'Colormap',colormap('parula'), 'DataAspectRatio', [1 1 1]};
            set(gca, ax_opt{:});
            set(gca,'Visible','off')
            colorbar
            
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