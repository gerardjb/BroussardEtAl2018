function mask = maskUpdate(handles)
%Acts as the function callback for threshPick. Applies slider value to make
%and apply im2bw level, sets numeric edit display, and updates im2bw image
%in GUI.

global Img
global mask
global level

level = get(handles.Level,'Value');

mask = im2bw(Img,level);


set(handles.Display,'String',num2str(level));

axes(handles.axes1);
imshow(mask);