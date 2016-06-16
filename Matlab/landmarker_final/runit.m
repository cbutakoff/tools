%this 2D landmarker allows only to displace point, not add new!
%therefore it requires a predefined template which is specified by the 1st
%parameter, see example in AR_groups.groups. It is a matrix where every row
%corresponds to a separate polyline of the contour. If the whole contour 
%is given by a vector, Col 1 and 2 - represent first and last element for 
%that polyline, col 3 = 0 for open and =1 for closed polyline. Col 4 does
%not serve any purpose in this code.
%
%The second parameter to lm() is the line specification for the matlab plot
%command
%
%The shaqpe loading and saving are adapted to VTK XML polydata. Modify the
%functions to use different format.
%
%To start editing - choose image and shape from the interface (at least one
%shape must exist in the file), not necessarily with correct placement of
%landmarks.
%
%Update - reloads shapes and images from the current folder. By default
%only png images are shown
%
%Clicking on the shape file wil copy that filename into the "save filename"
%box, clicking on the image will generate the name for the shape.
%
%Clicking on the image will select the corresponding shape, clicking on the
%shape will not change the selected image. So different image-shape pairs
%can be used. At least as the initial contour.
%
%The shapes and images are assumed to be ordered in the same way when
%sorted alphabetically.
%
% Landmarking steps:
% 1) Put at least one VTP file in the folder with images. Keep the shape
% names the same  as the image names.
% 2) Run the script, choose folder, or click update if the current folder
% is used
% 3) Select image #1, the only vtp file will be selected automatically,
% click open and move the points around
% 4) When inished, do not close any window! Swich to the image selection
% window
% 5) Click save (or modify the filename in the save box and click save).
% The default name is image_name.vtp. After saving the lists should update
% automatically.
% 6) Select image #2, there is still only one shape, it will stay
% selected, open, edit, save.
% 7) Continue until bored
%
% May require a couple of tries before you get used to the interface. If
% something fails, just close and open. ;)


load AR_landmarks_subgroups
lm(AR_groups.groups,'y-x')