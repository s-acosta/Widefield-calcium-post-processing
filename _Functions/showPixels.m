% Iniatiation

W = importdata('WTR040_201205_WF.mat');
N = NeurotarExtractor('WTR040_201205');
XY = [N.behavior.X; N.behavior.Y];
XY(:, ~N.behavior.moving_time) = [];
W.moving_time = N.behavior.moving_time;
W.computeDFF;
W.maskSession;

% We get only the left hemisphere pixels of SSC
ssc_pixels = W.dff;
SSC = W.mask.SSC;
SSC(:,180:end) = false;  
SSC = repmat(SSC, [1 1 size(W.dff,3)]);
ssc_pixels(~SSC) = 0;
SSC = reshape(SSC(:,:,1),400, 400);

% We plot 5 random pixels to see their response
n_pixels = 5;
ppc_pxl_val = zeros(2,n_pixels);
ppc_trues = find(PPC);

figure,
tiledlayout(5,1)
for i = 1:5
    
    element = randi(numel(ppc_trues));
    [ppc_pxl_val(1,i), ppc_pxl_val(2,i)] = ind2sub([400, 400], ppc_trues(element));
    
    nexttile
    plot(squeeze(dff(ppc_pxl_val(1,i), ppc_pxl_val(2,i), 1:1500)))
    hold on
    plot(20*moving_time(1:1500), 'LineWidth', 2)
    hold off
    
end

% For the 5 random pixels, we compute de deconvolution

dff_denoised = zeros(n_pixels, 1500);
spikes = zeros(n_pixels, 1500);
tic
for i = 1:n_pixels
    pixel_trace = squeeze(ssc_pixels(ssc_pxl_val(1,i), ssc_pxl_val(2,i), 1:1500));
    [dff_denoised(i,:), spikes(i,:) , options] = deconvolveCa(pixel_trace);
end 
toc

%Show both set 

figure('Position', [41, 102, 1382, 851])
tiledlayout(5,1)
for i = 1:5
    
    pixel_trace = squeeze(ssc_pixels(ssc_pxl_val(1,i), ssc_pxl_val(2,i), 1:1500));
   
    nexttile
    
%     plot(pixel_trace, 'b')
%     hold on
    plot(dff_denoised(i,:), 'r')
    hold on
    plot(spikes(i,:), 'k')
    hold on
    plot(25*N.behavior.moving_time(1:1500), 'LineWidth', 2)
    hold off
    
end

% De-noise the whole signal
[ssc_pixels, ssc_recovery] = W.zipMat(ssc_pixels);
ssc_pixels_denoised = ssc_pixels;

for i = 1:size(ssc_pixels, 1)
    ssc_pixels_denoised(i,:) = deconvolveCa(ssc_pixels(i,:));
    subroutine_progressbar(i/size(ssc_pixels, 1));
end

%Spatial Info denoised vs regular
dff_denoised = zeros(size(W.dff,1), size(W.dff,2), 'single');

for i = 1:size(dff_denoised, 1)
    dff_denoised(i,:) = deconvolveCa(W.dff(i,:));
    subroutine_progressbar(i / size(W.dff,1));
end


mask_region = W.mask.ALL;
mask_region(~W.mask.Window) = 0;

figure,
subplot(1,2,1)
imagesc(mean(W.dff, 3), 'AlphaData', mask_region)
colorbar

subplot(1,2,2)
imagesc(mean(dff_denoised_mat,3), 'AlphaData', mask_region)
colorbar



SI_px_bn_raw = alloSpaceIndex(XY, W.dff , mean(W.dff,2), 8);
SI_px_bn_denoised = alloSpaceIndex(XY, dff_denoised , mean(dff_denoised,2), 8);

SI_px_raw = sum(SI_px_bn_raw, 2);
SI_px_denoised = sum(SI_px_bn_denoised,2);
SI_mat_raw = W.zipMat(SI_px_raw, W.recovery_vec);
SI_mat_denoised = W.zipMat(SI_px_denoised, W.recovery_vec);

%

figure,
subplot(1,2,1)
imagesc(SI_mat_raw)
caxis([0 prctile(SI_mat_raw(:), 99)])
colorbar

subplot(1,2,2)
imagesc(SI_mat_denoised)
caxis([0 prctile(SI_mat_denoised(:), 99)])
colorbar

%% Will LR work for thresholded pixels?

mask = W.mask.ALL;
mask(~W.mask.Window) = 0;
W.unzipDFF;

[traces, pixel_idx] = W.getTraces(W.dff, mask, 10);
W.showRegionTraces(W.dff, mask, traces, pixel_idx);

traces_deconvoluted = zeros(size(traces,1), size(traces,2));
for i = 1:size(traces,1)
    traces_deconvoluted(i,:) = squeeze(W.dff_deconvolved(pixel_idx(i,1), pixel_idx(i,2), :));
end

traces_deconvoluted = zeros(size(traces,1), size(traces,2));

for i = 1:10
    
    traces_deconvoluted(i,:) = lucric(traces(i,:)', 0.95, 2, 500);
end

%%
figure('Position', [200 200 800 400])
tiledlayout(10,2)
for i = 1:size(traces,1)
    nexttile(2*i - 1)
    plot(traces(i,:));
    nexttile(2*i)
    plot(traces_deconvoluted(i,:))
end










