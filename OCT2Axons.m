function varargout = OCT2Axons(varargin)
%USAGE:
%    [Axons_segmented] = OCT2Axons(OCTPath, AxonThr, Min_length)
%INPUTS:
%   - OCTPath - Input en face images at .tiff stack
%   - AxonThr - 0-100 value of thresholding an Axon after convulotion 
%   - Min_length - Minimum length of an Axon in pixels
%OUTPUTS:
%   - Axons_segmented - 3D Axons mask

tic
%% Load & Init paramaters
    Original_tiff_info = imfinfo(OCTpath);
    Original_tiff_size = size(Original_tiff_info,1);
    [Original] = uint8((zeros(Original_tiff_info(1).Height,Original_tiff_info(1).Width,Original_tiff_size)));
    [Seg_topBottom,Seg_leftRight,Seg_diagonal1,...
        Seg_diagonal2] = deal((zeros(Original_tiff_info(1).Height,Original_tiff_info(1).Width,Original_tiff_size)));
    
    for i=1:Original_tiff_size
        Original(:,:,i) = imread(OCTpath,i);
    end
fprintf('Loading run time: %d min and %.1f seconds\n', fix(toc/60),rem(toc,60))
tic

%%
F = [-1 -1 -1 1 1 1 -1 -1 -1]; %Hadamard inspired column vector
Seg_matrix = repmat(F,size(F,2),1);
Seg_matrix_diagonal1 = imrotate(Seg_matrix,45,'bilinear','crop');
Seg_matrix_topBottom = Seg_matrix-mean(Seg_matrix(:));
Seg_matrix_topBottom = Seg_matrix_topBottom./sqrt(sum(Seg_matrix_topBottom(:).^2));
Seg_matrix_leftRight = Seg_matrix_topBottom';
Seg_matrix_diagonal1 = Seg_matrix_diagonal1-mean(Seg_matrix_diagonal1(:));
Seg_matrix_diagonal1 = Seg_matrix_diagonal1./sqrt(sum(Seg_matrix_diagonal1(:).^2));
Seg_matrix_diagonal2 = imrotate(Seg_matrix,-45,'bilinear','crop');
Seg_matrix_diagonal2 = Seg_matrix_diagonal2-mean(Seg_matrix_diagonal2(:));
Seg_matrix_diagonal2 = Seg_matrix_diagonal2./sqrt(sum(Seg_matrix_diagonal2(:).^2));

%%
    for i=1:size(Original,3)
        Seg_topBottom_pre = conv2(Original(:,:,i),Seg_matrix_topBottom,'same');
        Seg_leftRight_pre = conv2(Original(:,:,i),Seg_matrix_leftRight,'same');
        Seg_diagonal1_pre = conv2(Original(:,:,i),Seg_matrix_diagonal1,'same');
        Seg_diagonal2_pre = conv2(Original(:,:,i),Seg_matrix_diagonal2,'same');
        Seg_topBottom(:,:,i) = Seg_topBottom_pre > prctile(Seg_topBottom_pre,AxonThr,'all');
        Seg_leftRight(:,:,i) = Seg_leftRight_pre > prctile(Seg_leftRight_pre,AxonThr,'all');
        Seg_diagonal1(:,:,i) = Seg_diagonal1_pre > prctile(Seg_diagonal1_pre,AxonThr,'all');
        Seg_diagonal2(:,:,i) = Seg_diagonal2_pre > prctile(Seg_diagonal2_pre,AxonThr,'all');
    end

Axons_segmented_bw = (false(size(Seg_topBottom ,1),size(Seg_topBottom ,2),size(Seg_topBottom,3)));
    for i=1:size(Seg_topBottom,3)
        Axons_segmented_pre = Seg_topBottom(:,:,i) | Seg_leftRight(:,:,i) | Seg_diagonal1(:,:,i) | Seg_diagonal2(:,:,i);
        Axons_segmented_bw(:,:,i) = bwareaopen(Axons_segmented_pre,Min_length,6);
    end
fprintf('Axons Segmentation run time: %d min and %.1f seconds\n', fix(toc/60),rem(toc,60))

%%
varargout{1} = Axons_segmented_bw;
end