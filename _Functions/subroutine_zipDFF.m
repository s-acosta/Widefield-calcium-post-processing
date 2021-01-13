function [DFF_new, recovery_vec] = subroutine_zipDFF(DFF_old, recovery_vec)
% The goal of this subroutine is to reduce the size of large matrices that
% have been masked with zeros. It removes the zeros so that it can be saved
% (along with the recovery vector). If the recovery vector is introduced,
% then it takes the matrix back to its normal shape. The function assumes
% the input matrices to be either [x_pixels y_pixels n_frames] or [n_pixels
% n_frames].
%
% Developed by Santi on 29/12/2020
% Last update 29/12/2020

if nargin == 1
    disp('Zipping mode')
    
    n_pixels = size(DFF_old, 1) * size(DFF_old,2);
    n_frames = size(DFF_old, 3);
    
    DFF_old = reshape(DFF_old, [n_pixels, n_frames]);
    recovery_vec = DFF_old(:,1) ~= 0;
    
    DFF_new = DFF_old;
    DFF_new(~recovery_vec, :) = [];
    
elseif nargin == 2
    disp('Unzipping mode')
    
    n_pixels = length(recovery_vec);
    n_frames = size(DFF_old, 2);
    
    DFF_new = zeros(n_pixels, n_frames, 'single');
    
    k = 0;
    for i = 1:n_pixels
        
        if recovery_vec(i)
            k = k + 1;
            DFF_new(i,:) = DFF_old(k,:);
        end
        
    end
    
    DFF_new = reshape(DFF_new, [sqrt(n_pixels) sqrt(n_pixels) n_frames]);
    
end

end