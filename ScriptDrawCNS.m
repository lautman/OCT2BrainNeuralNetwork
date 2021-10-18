function [circlePixels] = ScriptDrawCNS(Img,S_CNS_segmented_cell_pre,CellRadius,MinzLength)
%USAGE:
%   [circlePixels] = ScriptDrawCNS(Img,S_CNS_segmented_cell_pre,CellRadius)
%INPUTS:
%   - Img - (pre) single CNS mask en face image
%   - S_CNS_segmented_cell_pre - regionprops of Img
%   - CellRadius - CNS radius (RA)
%   - MinzLength - Minimum number of slices a CNS could spread
%OUTPUTS:
%   - circlePixels - (post) CNS mask en face image

imageSizeX = size(Img,1);
imageSizeY = size(Img,2);
[columnsInImage,rowsInImage] = meshgrid(1:imageSizeX, 1:imageSizeY);
circlePixels = uint8(zeros(size(Img,1),size(Img,2)));
    for i=1:size(S_CNS_segmented_cell_pre,1)
    centerX = round(S_CNS_segmented_cell_pre(i).Centroid(1));
    centerY = round(S_CNS_segmented_cell_pre(i).Centroid(2));
    radius = ceil(S_CNS_segmented_cell_pre(i).MaxFeretDiameter./2);
        if (radius >= MinzLength) % cell size cut off 
            radius = CellRadius;
            circlePixels = circlePixels + uint8((rowsInImage - centerY).^2 ...
                + (columnsInImage - centerX).^2 <= radius.^2).*255;
        end
    end
end