classdef MaskGUI < matlab.apps.AppBase
    
    % obj = MaskGUI() class creates a GUI that matches an already created
    % template with an image. This class does so by zooming and translating
    % the template matrix so that it is alineated with the image. The
    % overall goal is to mask the original image in different regions.
    %
    % MaskGUI Properties:
    %    template - Original template
    %    image - Original image
    %    move_order - Transformations applied to the template
    %    curr_template - Modified template
    %
    % MaskGUI Methods:
    %    doThis - Description of doThis
    %    doThat - Description of doThat
    
    % Properties that correspond to app components
    properties (Access = private)
        UIFigure            matlab.ui.Figure
        ImageAxes           matlab.ui.control.UIAxes
        LoadImageButton     matlab.ui.control.Button
        LoadTemplateButton  matlab.ui.control.Button
        UpButton            matlab.ui.control.Button
        LeftButton          matlab.ui.control.Button
        RightButton         matlab.ui.control.Button
        DownButton          matlab.ui.control.Button
        ZoomInButton        matlab.ui.control.Button
        ZoomOutButton       matlab.ui.control.Button
        DropDown            matlab.ui.control.DropDown
        WindowButton        matlab.ui.control.Button
        FinishButton        matlab.ui.control.Button
        
    end
    
    % Properties that correspond to the template and its image
    properties (Access = public)
        template
        curr_template
        window_mask
        image
    end
    
    properties (Access = private)
        move_order
        mask_order
        curr_image
    end
    
    % Loading and showing related methods
    methods (Access = private)
        
        function loadImage(app,image_file)
            
            try
                app.image = imread(image_file);
                
            catch ME
                % If problem reading image, display error message
                uialert(app.UIFigure, ME.message, 'Image Error');
                return;
            end
            
            try
                clearImage(app)
            catch
            end
            
        end
        
        function showImage(app)
            
            clearImage(app)
            maskImage(app);
            imagesc(app.ImageAxes, app.curr_image);
            
        end
        
        function clearImage(app)
            
            plots = get(app.ImageAxes,'Children');
            
            for i = 1:length(plots)
                
                if contains(class(plots(i)), 'Image')
                    
                    delete(app.ImageAxes.Children(i))
                    
                end
                
            end
            
        end
        
        function loadTemplate(app, template_file)
            
            try
                app.template = importdata(template_file);
                app.template = bwperim(app.template, 8);
                
            catch ME
                % If problem reading template, display error message
                uialert(app.UIFigure, ME.message, 'Template Error');
                return;
            end
            
        end
        
        function showTemplate(app)
            
            clearTemplate(app)
            moveTemplate(app)
            
            hold(app.ImageAxes, 'on')
            
            template_2show = app.curr_template;
            for i = 1:size(app.curr_template, 3)
                single_layer = template_2show(:, :, i);
                single_layer(~app.window_mask) = 0;
                template_2show(:, :, i) = single_layer;
            end
            
            template_colors = 0.85  * colormap(hsv(7));
            
            if isequal(app.mask_order, 0) || isempty(app.mask_order)
                
                for i = 1:size(app.curr_template, 3)
                    contour(app.ImageAxes, template_2show(:,:,i), [1 1],...
                        'Color',template_colors(i,:),'LineWidth', 2);
                end
                
            else
                
                contour(app.ImageAxes, template_2show(:,:, app.mask_order), ...
                    [1 1],...
                    'Color',template_colors(app.mask_order,:),...
                    'LineWidth', 2);
                
            end
            
            hold(app.ImageAxes, 'off')
            
        end
        
        function clearTemplate(app)
            
            plots = get(app.ImageAxes,'Children');
            true_plot = zeros(1,length(plots));
            
            for i = 1:length(true_plot)
                
                if contains(class(plots(i)), 'Contour')
                    
                    true_plot(i) = i;
                    
                else
                    true_plot(i) = [];
                    
                end
                
                
            end
            
            delete(app.ImageAxes.Children(true_plot))
            
        end
        
    end
    
    % Template moving methos
    methods (Access = private)
        
        function registerOrder(app, order)
            
            switch order
                
                case 'Right'
                    app.move_order(1) =  app.move_order(1) + 1;
                    
                case 'Left'
                    app.move_order(1) =  app.move_order(1) - 1;
                    
                case 'Up'
                    app.move_order(2) =  app.move_order(2) + 1;
                    
                case 'Down'
                    app.move_order(2) =  app.move_order(2) - 1;
                    
                case 'In'
                    app.move_order(3) = app.move_order(3) + 1;
                    
                case 'Out'
                    app.move_order(3) = app.move_order(3) - 1;
                    
            end
            
            
        end
        
        function moveTemplate(app)
            
            % moveTemplate  Transforms the template to match the image
            %   The goal of this function is to apply the already specified
            %   transformations of registerOrder to the original template.
            %
            %   Every time it starts from scratch, so that whatever the
            %   former transformations is, it won't affect the following (e.g. if
            %   the template has moved out of the boundaries in the
            %   previous transformations we are not losing that
            %   information).
            %
            %   The transformations are specified by the array move_order.
            %   move_order(1) specifices lateral translation, with negative
            %   being left and positive being right. Similarly,
            %   move_order(2) follows the same logic for up and down and
            %   move_order(3) does the appropiate for zoomin/zoomout, with
            %   zooming being made with resizing of the order of 0.025
            %
            %   In order to do so, frame is first created. Frame is a
            %   matrix 3 times the template. The goal of frame is not to
            %   lose information if multiple orders are applied. After all
            %   transformations are applied the template (center of matrix frame)
            %   is retrieved.
            %
            %   See also registerOder.
            
            app.curr_template = app.template;
            
            if nnz(app.move_order) == 0
                return
            end
            
            n_row = size(app.curr_template,1);
            n_column = size(app.curr_template, 2);
            n_layers = size(app.curr_template, 3);
            
            frame = zeros(3*n_row, 3*n_column, n_layers);
            frame_centre = [floor(3*n_row/2), floor(round(3*n_row/2))];
            
            % Zoom in/out
            new_template = imresize(app.curr_template, ...
                1 + app.move_order(3) * 0.025);
            
            is_odd = mod(size(new_template, [1 2]), 2) == 1;
            template_centre = floor(size(new_template, [1 2]) / 2);
            
            % Where to put the template inside the frame
            frame_start = frame_centre - template_centre;
            frame_end = frame_centre + template_centre;
            frame_vec = [frame_start(1)+1 : frame_end(1) + is_odd(1); ...
                frame_start(2)+1 : frame_end(2) + is_odd(2)];
            
            frame(frame_vec(1,:), frame_vec(2,:), :) = new_template;
            
            % Move right/left & up/down
            frame = circshift(frame, [-app.move_order(2) app.move_order(1) 0]);
            
            % Extract the template from the frame
            app.curr_template = frame(n_row+1 : 2*n_row, n_column+1 : 2*n_column, :);
            
        end
        
    end
    
    % Template masking methods
    methods (Access = private)
        
        function registerArea(app, area)
            
            switch area
                
                case 'All', app.mask_order = 0;
                    
                case 'SMC', app.mask_order = 1;
                    
                case 'SSC', app.mask_order = 2;
                    
                case 'PPC', app.mask_order = 3;
                    
                case 'RSC', app.mask_order = 4;
                    
                case 'VC', app.mask_order = 5;
                    
                case 'AC', app.mask_order = 6;
                    
                case 'TAC', app.mask_order = 7;
                    
                otherwise, app.mask_order = [];
                    
            end
            
        end
        
        function maskImage(app)
            
            app.curr_image = app.image;
            
            if ~isempty(app.window_mask)
                app.curr_image(~app.window_mask) = 0;
            end
            
            if  isempty(app.mask_order)
                return
            else
                mask = createTemplateMask(app, app.mask_order);
            end
            
            app.curr_image(~mask) = 0;
           
        end
        
        
        function mask = createTemplateMask(app, layer)
            
            if layer == 0
                mask = sum(app.curr_template, 3);
            else
                mask = app.curr_template(:, :, layer);
            end
            
            mask = bwperim(mask, 8);
            
            % Make top line white:
            col1 = find(mask(1, :), 1, 'first');
            col2 = find(mask(1, :), 1, 'last');
            mask(1, col1:col2) = true;
            
            % Make bottom line white:
            col1 = find(mask(end, :), 1, 'first');
            col2 = find(mask(end, :), 1, 'last');
            mask(end, col1:col2) = true;
            
            % Make right line white:
            row1 = find(mask(:, end), 1, 'first');
            row2 = find(mask(:, end), 1, 'last');
            mask(row1:row2, end) = true;
            
            % Make left line white:
            row1 = find(mask(:, 1), 1, 'first');
            row2 = find(mask(:, 1), 1, 'last');
            mask(row1:row2, 1) = true;
            
            mask = imfill(mask, 'holes');
            
        end
        
    end
    
    % Window masking methods
    methods (Access = private)
        
        function windowMask(app)

            [rows, columns, ~] = size(app.curr_image);

            clearTemplate(app)
            
            try
                c = drawcircle(app.ImageAxes, 'Center',c.Center,...
                    'Radius',c.Radius, ...
                    'Color', 'w', 'FaceAlpha', 0.4);
            catch
                c = drawcircle(app.ImageAxes, 'Color', 'w', 'FaceAlpha', 0.4);
                
            end
            
            MaskGUI.customWait(c);
            
            angles = linspace(0, 2*pi, 10000);
            window_radius = c.Radius;
            window_center = c.Center;
            
            x = cos(angles) * window_radius + window_center(1);
            y = sin(angles) * window_radius + window_center(2);
            app.window_mask = poly2mask(x, y, rows, columns);
            
            showImage(app)

        end
        
    end
    
    methods (Access = private, Static = true)
        
        function pos = customWait(hROI)
            
            % Listen for mouse clicks on the ROI
            l = addlistener(hROI,'ROIClicked',@MaskGUI.clickCallback);
            
            % Block program execution
            uiwait;
            
            % Remove listener
            delete(l);
            
            % Return the current position
            pos = hROI.Position;
            
        end
        
        function clickCallback(~,evt)
            
            if strcmp(evt.SelectionType,'double')
                uiresume;
            end
            
        end
        
        
    end
    
    % Callbacks that handle component events
    methods (Access = private)
        
        function startupFcn(app)
            
            % Mock image to start with
            if isempty(app.image)
                loadImage(app, 'WTR039_201116_refimg.tif');
            end
            
            % Mock template to start with
            if isempty(app.template)
                loadTemplate(app, 'ROIs.mat');
            end
            
            showImage(app);
            showTemplate(app)
            
        end
        
        function LoadImageButtonPushed(app, event)
            
            % Display uigetfile dialog
            filterspec = {'*.tif;*.png;*.gif','All Image Files'};
            [f, p] = uigetfile(filterspec);
            
            % Make sure user didn't cancel uigetfile dialog
            if (ischar(p))
                fname = [p f];
                loadImage(app, fname);
            end
        end
        
        function LoadTemplateButtonPushed(app, event)
            
            % Display uigetfile dialog
            filterspec = {'*.mat'};
            [f, p] = uigetfile(filterspec);
            
            % Make sure user didn't cancel uigetfile dialog
            if (ischar(p))
                fname = [p f];
                loadTemplate(app, fname);
            end
        end
        
        function LeftButtonPushed(app, event)
            
            registerOrder(app, 'Left')
            showImage(app)
            showTemplate(app)
            
        end
        
        function RightButtonPushed(app, event)
            
            registerOrder(app, 'Right')
            showImage(app)
            showTemplate(app)
            
        end
        
        function UpButtonPushed(app, event)
            
            registerOrder(app, 'Up')
            showImage(app)
            showTemplate(app)
            
        end
        
        function DownButtonPushed(app, event)
            
            registerOrder(app, 'Down')
            showImage(app)
            showTemplate(app)
            
        end
        
        function ZoomInButtonPushed(app, event)
            
            registerOrder(app, 'In')
            showImage(app)
            showTemplate(app)
            
        end
        
        function ZoomOutButtonPushed(app, event)
            
            registerOrder(app, 'Out');
            showImage(app)
            showTemplate(app)
            
        end
        
        function DropDownValueChanged(app, event)
            
            registerArea(app, app.DropDown.Value)
            showImage(app)
            showTemplate(app)
            
        end
        
        function WindowButtonPushed(app, event)
            
            windowMask(app)
            
        end
        
        function FinishButtonPushed(app, event)
            
            Output = struct;
            Output.Window = app.window_mask;
            Output.ALL = createTemplateMask(app, 0);
            Output.SMC = createTemplateMask(app, 1);
            Output.SSC = createTemplateMask(app, 2);
            Output.PPC = createTemplateMask(app, 3);
            Output.RSC = createTemplateMask(app, 4);
            Output.VC = createTemplateMask(app, 5);
            Output.AC = createTemplateMask(app, 6);
            Output.TAC = createTemplateMask(app, 7);
            
            assignin('base', 'mask', Output)
            
            delete(app)
            
        end
        
    end
    
    % Component initialization
    methods (Access = private)
        
        % Create UIFigure and components
        function createComponents(app)
            
            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.AutoResizeChildren = 'off';
            app.UIFigure.Position = [100 100 1000 750];
            app.UIFigure.Name = 'Template Matcher';
            app.UIFigure.Resize = 'off';
            app.UIFigure.Color = [0 0 0];
            
            % Create UIAxes
            app.ImageAxes = uiaxes(app.UIFigure);
            axis(app.ImageAxes, 'image')
            app.ImageAxes.Color = [0 0 0];
            app.ImageAxes.XTick = [];
            app.ImageAxes.XTickLabel = {'[ ]'};
            app.ImageAxes.YTick = [];
            app.ImageAxes.Position = [334 66 631 562];
            app.ImageAxes.XColor = [0 0 0];
            app.ImageAxes.YColor = [0 0 0];
            app.ImageAxes.Colormap = colormap('gray');
            app.ImageAxes.BackgroundColor = [0 0 0];
            
            % Create LoadImageButton
            app.LoadImageButton = uibutton(app.UIFigure, 'push');
            app.LoadImageButton.Position = [38 572 115 22];
            app.LoadImageButton.Text = 'Load Image';
            app.LoadImageButton.ButtonPushedFcn = createCallbackFcn(app, @LoadImageButtonPushed, true);
            
            % Create LoadTemplateButton
            app.LoadTemplateButton = uibutton(app.UIFigure, 'push');
            app.LoadTemplateButton.Position = [197 572 115 22];
            app.LoadTemplateButton.Text = 'Load Template';
            app.LoadTemplateButton.ButtonPushedFcn = createCallbackFcn(app, @LoadTemplateButtonPushed, true);
            
            % Create LeftButton
            app.LeftButton = uibutton(app.UIFigure, 'push');
            app.LeftButton.Position = [53 434 100 25];
            app.LeftButton.Text = 'Left';
            app.LeftButton.ButtonPushedFcn = createCallbackFcn(app, @LeftButtonPushed, true);
            
            % Create RightButton
            app.RightButton = uibutton(app.UIFigure, 'push');
            app.RightButton.Position = [197 434 100 25];
            app.RightButton.Text = 'Right';
            app.RightButton.ButtonPushedFcn = createCallbackFcn(app, @RightButtonPushed, true);
            
            % Create UpButton
            app.UpButton = uibutton(app.UIFigure, 'push');
            app.UpButton.Position = [129 475 100 25];
            app.UpButton.Text = 'Up';
            app.UpButton.ButtonPushedFcn = createCallbackFcn(app, @UpButtonPushed, true);
            
            % Create DownButton
            app.DownButton = uibutton(app.UIFigure, 'push');
            app.DownButton.Position = [129 386 100 25];
            app.DownButton.Text = 'Down';
            app.DownButton.ButtonPushedFcn = createCallbackFcn(app, @DownButtonPushed, true);
            
            % Create ZoomInButton
            app.ZoomInButton = uibutton(app.UIFigure, 'push');
            app.ZoomInButton.Position = [53 307 100 22];
            app.ZoomInButton.Text = 'Zoom In';
            app.ZoomInButton.ButtonPushedFcn = createCallbackFcn(app, @ZoomInButtonPushed, true);
            
            % Create ZoomOutButton
            app.ZoomOutButton = uibutton(app.UIFigure, 'push');
            app.ZoomOutButton.Position = [197 307 100 22];
            app.ZoomOutButton.Text = 'Zoom Out';
            app.ZoomOutButton.ButtonPushedFcn = createCallbackFcn(app, @ZoomOutButtonPushed, true);
            
            % Create DropDown
            app.DropDown = uidropdown(app.UIFigure);
            app.DropDown.Items = {'No mask','All', 'SMC', 'SSC', 'PPC', 'RSC', 'VC', 'AC', 'TAC'};
            app.DropDown.FontSize = 14;
            app.DropDown.Position = [120 230 109 28];
            app.DropDown.Value = 'No mask';
            app.DropDown.ValueChangedFcn = createCallbackFcn(app, @DropDownValueChanged, true);
            
            % Create Button
            app.WindowButton = uibutton(app.UIFigure, 'push');
            app.WindowButton.Position = [117 155 124 33];
            app.WindowButton.Text = 'Window Mask';
            app.WindowButton.ButtonPushedFcn = createCallbackFcn(app, @WindowButtonPushed, true);
            
            % Create FinishButton
            app.FinishButton = uibutton(app.UIFigure, 'push');
            app.FinishButton.Position = [102 77 153 40];
            app.FinishButton.Text = 'Finish';
            app.FinishButton.ButtonPushedFcn = createCallbackFcn(app, @FinishButtonPushed, true);
            
            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
            
        end
    end
    
    % App creation and deletion
    methods (Access = public)
        
        % Construct app
        function app = MaskGUI(image_file, template_file)
            
            app.move_order = zeros(1,3);
            app.mask_order = [];
            
            % Create UIFigure and components
            createComponents(app)
            
            % Register the app with App Designer
            registerApp(app, app.UIFigure)
            
            try
                loadImage(app, image_file)
                loadTemplate(app, template_file)
            catch
            end
            
            % Execute the startup function
            runStartupFcn(app, @startupFcn)
            
            if nargout == 0
                clear app
            end
        end
        
        % Code that executes before app deletion
        function delete(app)
               
            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end