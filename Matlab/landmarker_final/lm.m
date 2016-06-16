function varargout = lm(varargin)
% !!!first parameter - landmark groups, then whatever you want to pass to 
% !!!shape plotting function
%
% LM M-file for lm.fig
%      LM, by itself, creates a new LM or raises the existing
%      singleton*.
%
%      H = LM returns the handle to a new LM or the handle to
%      the existing singleton*.
%
%      LM('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in LM.M with the given input arguments.
%
%      LM('Property','Value',...) creates a new LM or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before lm_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to lm_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Copyright 2002-2003 The MathWorks, Inc.

% Edit the above text to modify the response to help lm

% Last Modified by GUIDE v2.5 10-Apr-2006 17:33:09

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @lm_OpeningFcn, ...
                   'gui_OutputFcn',  @lm_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT




% --- Executes just before lm is made visible.
function lm_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to lm (see VARARGIN)

% Choose default command line output for lm
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes lm wait for user response (see UIRESUME)
% uiwait(handles.figure1);

global the_image;

global fig;
fig = figure;
hold on;

global shift_pressed;
global control_pressed;
global alt_pressed;
global pan_pressed;

global undo_history; %one level currently

shift_pressed = false;
control_pressed = false;
alt_pressed = false;
pan_pressed = false;

%the_point = [100,100,200,200,100,200];
%the_point1 = [300,300,340,340,360,360];
%plot( the_point(1:2:end), the_point(2:2:end), '-yo' );
%plot( the_point1(1:2:end), the_point1(2:2:end), '-ro' );


set( fig, 'WindowButtonDownFcn', @MLBDown);
set( fig, 'WindowButtonUpFcn', @MLBUp);
set( fig, 'WindowButtonMotionFcn', @MLBMove);
set( fig, 'DoubleBuffer', 'on' );
set( fig, 'KeyPressFcn', @FigureKeyPress );
set( fig, 'KeyReleaseFcn', @FigureKeyRelease );


global move_started; %mouse clicked
global move_shape; %clicked in the vicinity of a landmark
move_started = false; 
move_shape = false;

global the_line;
global mouse_pointer;
mouse_pointer = get( gcf, 'Pointer' );


global landmark_groups;
landmark_groups = varargin{1};

global plot_data;
plot_data = varargin(2:end);

last_path = '';
if exist('lm.mat')
    load lm last_path
end;
set( handles.edtPath, 'String', last_path );

if exist(last_path)
    btnUpdate_Callback(handles.btnUpdate, eventdata, handles);
end;






% --- Outputs from this function are returned to the command line.
function varargout = lm_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;





function MLBDown(src,eventdata)
global move_started;
global move_shape;
global pan_pressed;
move_started = true;

%find closest point
click_pt = get( gca, 'CurrentPoint' );
click_pt = [click_pt(1,1); click_pt(1,2)];

lines = findobj( 'type', 'line' );

global the_line;
min_dist = 1.e+100;
for i = 1:length(lines)
    line1.handle = lines(i);
    line1.Data = [get(line1.handle,'XData'); get(line1.handle,'YData')];
    dist = sum( (repmat(click_pt,1, size(line1.Data,2)) - line1.Data).^2 );
    [dst, line1.pt_index] = min(dist);
    
    if dst<min_dist
        the_line = line1;
        min_dist = dst;
    end;
    
end;

if min_dist>100 %squared
    move_shape = false;
else 
    move_shape = true;
end;

if ( ~pan_pressed ) && move_shape
    SaveUndo();
end;


global mouse_pointer;
if move_shape
    set( gcf, 'PointerShapeCData', ones(16, 16)*NaN);
    set( gcf, 'Pointer', 'custom' );
end;

%save the clicked point
global clicked_point;
pt = get( gca, 'CurrentPoint' );
clicked_point = [pt(1,1); pt(1,2)];




function MLBUp(src,eventdata)
global move_started;
global move_shape;

global mouse_pointer;
if move_shape
    set( gcf, 'Pointer', mouse_pointer );
end;

move_started = false;
move_shape = false;


function MLBMove(src,eventdata)
global move_started;
global shift_pressed;
global control_pressed;
global alt_pressed;
global pan_pressed;
global move_shape;
global clicked_point;


if move_started

    pt = get( gca, 'CurrentPoint' );
    pt = [pt(1,1); pt(1,2)];
   
    if (~pan_pressed) && move_shape 


        global the_line;

        %if it is the first or the last, then if the ends are joined - move
        %them together.
        move2 = false;
        ind2 = 0;

        last_ind = size(the_line.Data,2);

        if( the_line.pt_index==1 )
            if( ~any( the_line.Data(:,1) ~= the_line.Data(:,last_ind) )  )
                move2 = true;
                ind2 = last_ind;
            end;
        elseif ( the_line.pt_index == last_ind )
            if( ~any( the_line.Data(:,1) ~= the_line.Data(:,last_ind) ) )
                move2 = true;
                ind2 = 1;
            end;
        end;


        old_pt = the_line.Data(:,the_line.pt_index);
        diplacement = pt-old_pt;

        if shift_pressed
            the_line.Data = the_line.Data + repmat(diplacement,1,size(the_line.Data,2));
        else
            the_line.Data(:,the_line.pt_index) = pt;

            if move2
                the_line.Data(:,ind2) = pt;
            end;
        end

        if control_pressed %this ignores the shift
            lines = findobj( 'type', 'line' );

            for i = 1:length(lines)
                line1.handle = lines(i);
                line1.Data = [get(line1.handle,'XData'); get(line1.handle,'YData')];
                line1.Data = line1.Data + repmat(diplacement,1,size(line1.Data,2));
                set( line1.handle, 'XData', line1.Data(1,:), 'YData', line1.Data(2,:) );

                %update also stored line
                if line1.handle == the_line.handle
                    line1.pt_index = the_line.pt_index; %this is my field
                    the_line = line1;
                end;
            end;

        else
            set( the_line.handle, 'XData', the_line.Data(1,:), 'YData', the_line.Data(2,:) );
        end;
    elseif pan_pressed %% do panning
        displ = clicked_point-pt;
        lims = [xlim; ylim];

        newlim = lims+repmat(displ,1,2);
        xlim(newlim(1,:));
        ylim(newlim(2,:));
        
    end;
end;


% --- Executes on selection change in lbImages.
function lbImages_Callback(hObject, eventdata, handles)
% hObject    handle to lbImages (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns lbImages contents as cell array
%        contents{get(hObject,'Value')} returns selected item from lbImages
sel = get( hObject, 'Value' );
strings = get( hObject, 'String' );
imagefile = strings{sel};

dot_indices = find(imagefile=='.');
last_dot = dot_indices(end);
imagefile = [imagefile( 1:last_dot-1 ),'.vtp'];

set( handles.edtShapeName, 'String', imagefile );

nshapes = length( get( handles.lbShapes, 'String' ) );
if sel<=nshapes
    set( handles.lbShapes, 'Value', sel );
end;


% --- Executes during object creation, after setting all properties.
function lbImages_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lbImages (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in lbShapes.
function lbShapes_Callback(hObject, eventdata, handles)
% hObject    handle to lbShapes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns lbShapes contents as cell array
%        contents{get(hObject,'Value')} returns selected item from lbShapes
sel = get( hObject, 'Value' );
strings = get( hObject, 'String' );
imagefile = strings{sel};

set( handles.edtShapeName, 'String', imagefile );






% --- Executes during object creation, after setting all properties.
function lbShapes_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lbShapes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btnImageOpen.
function btnImageOpen_Callback(hObject, eventdata, handles)
% hObject    handle to btnImageOpen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
image_ind = get(handles.lbImages, 'Value' );
images = get( handles.lbImages, 'String' );
image_file = images{image_ind};

shape_ind = get(handles.lbShapes, 'Value' );
shapes = get( handles.lbShapes, 'String' );
shape_file = shapes{shape_ind};

path = get(handles.edtPath,'String');

global fig;
figure( fig );
clf;
the_image = imread( [path, image_file] );
imshow( the_image ); hold on;
Maximize( gcf );

global landmark_groups;
shape = LoadShape( [path, shape_file] );

global plot_data;
PlotShape( shape', landmark_groups, plot_data{:} );

SaveUndo();



function SaveUndo()

global undo_history;
undo_history = {};
lines = findobj( 'type', 'line' );
undo_history{length(lines)}=[];
for i = 1:length(lines)
    undo_history{i} = [get(lines(i),'XData'); get(lines(i),'YData')];
end;


function LoadUndo()

global undo_history;
lines = findobj( 'type', 'line' );
for i = 1:length(lines)
    line = undo_history{i};
    set(lines(i),'XData',line(1,:)); 
    set(lines(i),'YData',line(2,:));
end;










function edtPath_Callback(hObject, eventdata, handles)
% hObject    handle to edtPath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edtPath as text
%        str2double(get(hObject,'String')) returns contents of edtPath as a double


% --- Executes during object creation, after setting all properties.
function edtPath_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edtPath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btnUpdate.
function btnUpdate_Callback(hObject, eventdata, handles)
% hObject    handle to btnUpdate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

path = get( handles.edtPath, 'String' );
if path( length(path) ) ~= '/'
    path = [path, '/'];
    set( handles.edtPath, 'String', path );
end;

%get images
struct = [dir( [path,'*.tif'] );dir( [path,'*.jpg'] );dir( [path,'*.png'] );dir( [path,'*.jpeg'] )];

images = {};
for i=1:length( struct )
    images{i} = struct(i).name;
end;
images = sort(images);


%get shapes
struct = dir( [path,'*.vtp'] );

shapes = {};
for i=1:length( struct )
    shapes{i} = struct(i).name;
end;
shapes = sort(shapes);


set( handles.lbImages, 'String', images );
set( handles.lbShapes, 'String', shapes );





% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over btnUpdate.
function btnUpdate_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to btnUpdate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)




% --- Executes on key press over lbShapes with no controls selected.
function lbShapes_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to lbShapes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end







% --- Executes on button press in btnSaveShape.
function btnSaveShape_Callback(hObject, eventdata, handles)
% hObject    handle to btnSaveShape (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
lines = findobj( 'type', 'line' );

shape = [];

for i = length(lines):-1:1
    line1 = lines(i);
    shape1 = [get(line1,'XData'); get(line1,'YData')];

    last_ind = size( shape1,2 );
    
    %remove last point if it is same as the first
    %contour is closed
    remove_last = false;
    
    if last_ind>1
        if( ~any( shape1(:,1) ~= shape1(:,last_ind) ) )
            remove_last = true;
        end;
    end;
        
    if remove_last 
        shape = [shape, shape1(:,1:(last_ind-1) )];
    else
        shape = [shape, shape1];
    end;
end;

shape = reshape( shape, 1, numel(shape) );

path = get( handles.edtPath, 'String' );
shapename = get( handles.edtShapeName, 'String' );

overwrite = true;
if exist([path, shapename])
    k = menu(['File ', shapename, ' exists. Overwrite?'],'No','Yes');
    overwrite = k==2;
end;
    
if overwrite
    SaveShape( [path, shapename], shape );
end;

btnUpdate_Callback( handles.btnUpdate, eventdata, handles);


function edtShapeName_Callback(hObject, eventdata, handles)
% hObject    handle to edtShapeName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edtShapeName as text
%        str2double(get(hObject,'String')) returns contents of edtShapeName as a double


% --- Executes during object creation, after setting all properties.
function edtShapeName_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edtShapeName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in btnBrowse.
function btnBrowse_Callback(hObject, eventdata, handles)
% hObject    handle to btnBrowse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
folder = uigetdir( get( handles.edtPath,  'String' ));
if folder
    set( handles.edtPath, 'String', folder );
    set( handles.lbShapes, 'Value', 1 );    
    set( handles.lbImages, 'Value', 1 );    
    btnUpdate_Callback( handles.btnUpdate, eventdata, handles);
end;




function FigureKeyRelease(src,evnt)
global shift_pressed;
global control_pressed;
global alt_pressed;
global pan_pressed;


control_pressed = false;
alt_pressed = false;
shift_pressed = false;

for x=1:length(evnt.Modifier)
    switch(evnt.Modifier{x})
        case 'control'
            control_pressed = true;
        case 'alt'
            alt_pressed = true;
        case 'shift'
            shift_pressed = true;
    end
end

if strcmp( evnt.Key, 'q' )
    pan_pressed = false;
end





function FigureKeyPress(src,evnt)

lines = findobj( 'type', 'line' );

global shift_pressed;
global control_pressed;
global alt_pressed;
global pan_pressed;
global mouse_pointer;



for x=1:length(evnt.Modifier)
    switch(evnt.Modifier{x})
        case 'control'
            control_pressed = true;
        case 'alt'
            alt_pressed = true;
        case 'shift'
            shift_pressed = true;
    end
end


xystep = 10;
if shift_pressed 
    xystep=1;
end;


if strcmp( evnt.Key, 'rightarrow' )
    for i = 1:length(lines)
        line1 = lines(i);
        data = get(line1,'XData'); 
        data = data+xystep;
        set( line1, 'XData', data );
    end;

elseif strcmp( evnt.Key, 'leftarrow' )
    for i = 1:length(lines)
        line1 = lines(i);
        data = get(line1,'XData'); 
        data = data-xystep;
        set( line1, 'XData', data );
    end;

elseif strcmp( evnt.Key, 'uparrow' )
    for i = 1:length(lines)
        line1 = lines(i);
        data = get(line1,'YData'); 
        data = data-xystep;
        set( line1, 'YData', data );
    end;

elseif strcmp( evnt.Key, 'add' )
    for i = 1:length(lines)
        line1 = lines(i);
        data = get(line1,'YData'); 
        data = data*1.05;
        set( line1, 'YData', data );
        data = get(line1,'XData'); 
        data = data*1.05;
        set( line1, 'XData', data );
    end;
    
elseif strcmp( evnt.Key, 'subtract' )
    for i = 1:length(lines)
        line1 = lines(i);
        data = get(line1,'YData'); 
        data = data/1.05;
        set( line1, 'YData', data );
        data = get(line1,'XData'); 
        data = data/1.05;
        set( line1, 'XData', data );
    end;

    
elseif strcmp( evnt.Key, 'downarrow' )
    for i = 1:length(lines)
        line1 = lines(i);
        data = get(line1,'YData'); 
        data = data+xystep;
        set( line1, 'YData', data );
    end;
    
elseif strcmp( evnt.Key, 'y' ) | strcmp( evnt.Key, 'm' ) |...
        strcmp( evnt.Key, 'c' ) | strcmp( evnt.Key, 'r' ) |...
        strcmp( evnt.Key, 'g' ) | strcmp( evnt.Key, 'b' ) |...
        strcmp( evnt.Key, 'w' ) | strcmp( evnt.Key, 'k' ) 
    
    for i = 1:length(lines)
        line1 = lines(i);
        set( line1, 'Color', evnt.Key );
    end;

    
elseif strcmp( evnt.Key, 'p' )
    set( gcf, 'Pointer', mouse_pointer );

    
% elseif strcmp( evnt.Key, 'f' )
%     lines = findobj( 'type', 'line' );
% 
%     centr = [0;0];
%     npts = 0;
%     for i = 1:length(lines)
%         line1.handle = lines(i);
%         line1.Data = [get(line1.handle,'XData'); get(line1.handle,'YData')];
%       
%         centr = centr + sum( line1.Data, 2 );
%         npts = npts+size( line1.Data, 2 );
%     end;
%     centr = centr/npts;
%     
%     for i = 1:length(lines)
%         line1.handle = lines(i);
%         line1.Data = [get(line1.handle,'XData'); get(line1.handle,'YData')];
%       
%         line1.Data(1,:) = 2*centr(1)-line1.Data(1,:);
%         
%         set( line1.handle, 'XData', line1.Data(1,:), 'YData', line1.Data(2,:) );
%     end;

    
    

elseif strcmp( evnt.Key, 'q' )
    pan_pressed = true;
    
elseif strcmp( evnt.Key, 'u' )
    LoadUndo();
    
elseif strcmp( evnt.Key, 'a' )
    lims = [xlim; ylim];
    center = (lims(:,1)+lims(:,2))/2;
    extent = lims(:,2)-center;
    newlim = [center-extent/2,center+extent/2];
    xlim(newlim(1,:));
    ylim(newlim(2,:));
    
elseif strcmp( evnt.Key, 'z' )
    lims = [xlim; ylim];
    center = (lims(:,1)+lims(:,2))/2;
    extent = lims(:,2)-center;
    newlim = [center-extent*2,center+extent*2];
    xlim(newlim(1,:));
    ylim(newlim(2,:));

    
elseif strcmp( evnt.Key, 'v' )
    %get the shape
    lines = findobj( 'type', 'line' );

    shape = [];

    for i = length(lines):-1:1
        line1 = lines(i);
        shape1 = [get(line1,'XData'); get(line1,'YData')];

        last_ind = size( shape1,2 );

        %remove last point if it is same as the first
        %contour is closed
        remove_last = false;
        if( ~any( shape1(:,1) ~= shape1(:,last_ind) ) )
            remove_last = true;
        end;

        if remove_last 
            shape = [shape, shape1(:,1:(last_ind-1) )];
        else
            shape = [shape, shape1];
        end;
    end;

    shape = reshape( shape, 1, numel(shape) );
    
    load('AVCAR_Estim_b1_LR5.mat');
    num = AVCAR_Estim_b1_LR5(shape);
    
    strings = {...
        -0.34 -0.26 'L3';
        -0.24 -0.16 'L2';
        -0.14 -0.06 'L1';
        -0.04 0.04 'FR';
        0.06 0.14  'R1';
        0.16 0.24  'R2';
        0.26 0.34  'R3'};
    
    
    str = '';
    for kkk = 1:size(strings,1)
       if num>=strings{kkk,1} & num<=strings{kkk,2}
           str = strings{kkk,3};
           break;
       end;
    end;
        
    msgbox( [num2str(num),'; ', str] );
end;


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure

last_path = get( handles.edtPath, 'String' );
save lm last_path;


global fig;
valid=any(fig==findobj);
if valid 
    close(fig);
end;
delete(hObject);



% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over btnSaveShape.
function btnSaveShape_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to btnSaveShape (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


