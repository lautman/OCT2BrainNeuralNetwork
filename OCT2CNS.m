function varargout = OCT2CNS(varargin)
%USAGE:
%    [Seg_cell3, I_OCT] = OCT2CNS(OCTPath, MinCell, CellThr, [PixelX, PixelY, PixelZ],
%                                      MinzLength, RA);
%INPUTS:
%   - OCTPath - Input en face images at .tiff stack
%   - MinCell - Area of minimum CNS 
%   - CellThr - 0-100 value of thresholding a cell after convulotion 
%   - [PixelX, PixelY, PixelZ] - size of pixels in um
%   - MinzLength - Minimum number of slices a CNS could spread
%   - RA - CNS radius
%OUTPUTS:
%   - Seg_cell3 - 3D CNS mask
%   - I_OCT - Original Image 

tic
%% Input Checks & Init parameters
if (iscell(varargin{1}))
    %the first varible contains a cell with the rest of the varibles, open it
    varargin = varargin{1};
end 

OCTPath = varargin{1}; 
MinCell = varargin{2};
CellThr = varargin{3};
[PixelX, PixelY, PixelZ] = varargin{4};
MinzLength = varargin{5};
RA = varargin{6};

%% CNS Segmentation

% Filter Build
filt_init = zeros(RA * 4);
filt_init(round(RA*2),round(RA*2)) = 1;
Radius_circle1 = round(bwdist(filt_init));
Circle1 = ~(Radius_circle1 >= (round(RA*2)-1));
Circle2=(Radius_circle1<=(round(RA))).*-1;
Circle3=Circle2*2+Circle1;
filt_pre_final = Circle3-mean(Circle3(:));
filt_final = filt_pre_final./sqrt(sum(filt_pre_final(:).^2));

filt_mean = Circle2-mean(Circle2(:));
MinVolume = abs(sum(sum(Circle2)));

% Load OCT
Original_tiff_info = imfinfo(OCTPath);
Original_tiff_size = size(Original_tiff_info,1);
I_OCT = uint8((zeros(Original_tiff_info(1).Height,Original_tiff_info(1).Width,Original_tiff_size)));
I_OCT_Conv = zeros(Original_tiff_info(1).Height,Original_tiff_info(1).Width,Original_tiff_size);
Seg_cell1 = zeros(Original_tiff_info(1).Height,Original_tiff_info(1).Width,Original_tiff_size);
Seg_cell2 = zeros(Original_tiff_info(1).Height,Original_tiff_info(1).Width,Original_tiff_size);

for i=1:Original_tiff_size
    I_OCT(:,:,i) = imread(OCTPath,i);
end

fprintf('Loading run time: %d min and %.1f seconds\n', fix(toc/60),rem(toc,60))
tic

% Filter convolution
for i=1:size(I_OCT,3)
    I_OCT_Conv(:,:,i) = conv2(I_OCT(:,:,i),filt_final,'same');
    Seg_cell1(:,:,i) = I_OCT_Conv(:,:,i) > prctile(I_OCT_Conv(:,:,i),CellThr,'all');
    temp1 = bwareaopen(Seg_cell1(:,:,i),3);
    temp2 = regionprops(temp1,'Centroid','MaxFeretProperties');
    [circlePixels] = ScriptDrawCellBodies(temp1,temp2,RA,MinzLength);
    Seg_cell2(:,:,i) = bwareaopen(circlePixels,MinCell);
end

%% 3D Morphological Thresholding
MaxzLength = round(RA*2*PixelZ);
MinVolume = pi*(RA^2) * MinzLength  * 0.8;
MaxVolume = pi*(RA^2) * MaxzLength ; 
MinDiam = RA;
MaxDiam = RA*4;
Seg_cell2_Crop_CC = bwconncomp(Seg_cell2,6);
L = labelmatrix(Seg_cell2_Crop_CC);
S = regionprops3(Seg_cell2_Crop_CC,'VoxelList','Volume','image');
zLength = zeros(size(S,1),1);
zVolume = zeros(size(S,1),1);
zMaxFeretDiameter = zeros(size(S,1),1);
for j = 1:size(S,1)
   zLength(j) = round((S.VoxelList{j}(end,3)-S.VoxelList{j}(1,3)+1)*PixelZ);
   zVolume(j) = round(S.Volume(j)*PixelZ);
   temp = S.Image{j};
   tempMax = max(temp, [], 3);
   s = regionprops(tempMax,'MaxFeretProperties');
   zMaxFeretDiameter(j) = s.MaxFeretDiameter;
end

Seg_cell3 = ismember(L, find((zLength' >= MinzLength & zLength' <= MaxzLength &...
    zVolume' >= MinVolume & zVolume' <= MaxVolume &...
    zMaxFeretDiameter' >= MinDiam & zMaxFeretDiameter' <= MaxDiam))); % assumption cell spread through max 15um slice

fprintf('CNS Segmentation run time: %d min and %.1f seconds\n', fix(toc/60),rem(toc,60))
tic
    
%% Save varargout
varargout{1} = Seg_cell3;
varargout{2} = I_OCT;

end