% jimagesc New Version.
% Simple image analysis tool
%
% Usage: jimagesc(I, scale)
% 
% I     = image data, can be 3-dimensional matrix or cells. Imaginary data also available
% scale = scale factor, 1x2, can be omitted

 
function hhh=jimagesc_new(Ioriginal, scale)

hh=figure; 
set(hh, 'Toolbar','figure','visible','off','units', 'normalized', 'Position', [0.25, 0.25, 0.35, 0.5]);

handles=findall(hh,'type','uitoggletool'); 
index=logical(strcmp(get(handles,'Tag'),'Exploration.ZoomOut')+strcmp(get(handles,'Tag'),'Exploration.ZoomIn')); 
zoom_handles=handles(index); 

axis off;
AX='image';
CM='gray';
scaleLock=0;
zoomLock=0;
frameN=1;
sno = size(Ioriginal, 3);
FixedScale={};

% When 3-dimensional input    
if size(Ioriginal,3)>1
    MaxFrame=size(Ioriginal,3);
    type_text='matrix';
    Ibuffer=Ioriginal(:,:,frameN);
    I=Ibuffer;
    
% Cell type or 2-dimensional input        
elseif iscell(Ioriginal) ==1 && ( size(Ioriginal,2) >= 2 )
    MaxFrame=size(Ioriginal,2);
    type_text='cell';
    Ibuffer=Ioriginal{frameN};
    I=Ibuffer;
else
    Ibuffer=Ioriginal;
    I=Ioriginal;
end

%% For 3-Dimensional input, create necessary UI buttons.
if exist('type_text','var') 
    btn_prev = uicontrol('style','pushbutton','String','Prev Image','Callback', {@prev_pressed}, 'units','normalized','Position',[0.3 0.005 0.1 0.05] );
    btn_next = uicontrol('style','pushbutton','String','Next Image','Callback', {@next_pressed}, 'units','normalized','Position',[0.4 0.005 0.1 0.05] );
    scLock    = uicontrol('style','checkbox','String','Scale Lock','Callback', {@sc_Lock}, 'units','normalized','Position',[0.73 0.92 0.14 0.03], 'Fontsize', 10);
    zLock    = uicontrol('style','checkbox','String','Zoom Lock','Callback', {@zoom_Lock}, 'units','normalized','Position',[0.865 0.92 0.14 0.03], 'Fontsize', 10);
    frame_number = uicontrol('style','text','string',sprintf('Frame #: 1 / %d', MaxFrame), 'units','normalized','Position',[0.1 0.01 0.2 0.04], 'Fontsize', 10);
    sld = uicontrol('style', 'slider', 'SliderStep',[1/(sno-1) 10/(sno-1)], 'units', 'normalized', 'Position', [0.52, 0.01, 0.45, 0.04], 'Callback', @frame_shift);
    addlistener(sld,'ContinuousValueChange',@frame_shift);
end
%% In case of complex data input
ri_tag = uicontrol('style','popupmenu','String',{'absolute','phase','real','imaginary'},'Callback', {@ri_change}, 'units','normalized','Position',[0.3 0.86 0.1 0.1] );
set(ri_tag,'Visible','Off');
if isreal(I) ~= 1
    set(ri_tag,'Visible','On');
    I=abs(Ibuffer);
end
%% Automatically set scale if scale parameter is not assigned
if nargin==1
    scale=setscale(I);
end

%% Show the UI
h=axes('units','normalized','Position',[0.1 0.1 0.75 0.8]);
guard= abs(scale(2) - scale(1))/100;

% Scale Bar
slH=uicontrol('style','slider','units','normalized','Position',[0.87 0.1 0.03 0.8],'Max',scale(2), 'Min', scale(1), 'Value', scale(2), 'Callback', {@high_change});
slL=uicontrol('style','slider','units','normalized','Position',[0.93 0.1 0.03 0.8],'Max',scale(2), 'Min', scale(1), 'Value', scale(1), 'Callback', {@low_change} );

% Axis, Cmap Bar
axcontrol = uicontrol('style','popupmenu','String',{'image','tight'},'Callback', {@ratio_change}, 'units','normalized','Position',[0.1 0.86 0.1 0.1] );
c_map     = uicontrol('style','popupmenu','String',{'gray','jet','hot','parula'},'Callback', {@cmap_change}, 'units','normalized','Position',[0.2 0.86 0.1 0.1] );

% If the image has imaginary number, show the absolute valued image, first.
% Set ScrollWheel Function
set(hh,'windowscrollWheelFcn', @showImage);
imagesc(I, scale ); axis(AX); colormap(CM); colorbar

set(hh,'visible','on');

%% Define Callback Functions
    function high_change(source, eventdata)

        tmpH = get(slH,'Value');

        if scaleLock==0

            if tmpH > scale(1)
                scale(2)=tmpH;
%                 redraw;
            else       
                scale(2)= scale(1) + guard;
                set(slH,'Value',scale(2));
%                 redraw;
            end
            if zoomLock==1
                FixedScale = get(gca, {'xlim', 'ylim'});
            end
            redraw;
        end
    end

    function low_change(source, eventdata)
        
        tmpL = get(slL,'Value');

        if scaleLock==0
            
            if tmpL < scale(2)
                scale(1) = tmpL;
            else
                scale(1) = scale(2) - guard;
                set(slL,'Value',scale(1) );
            end
            if zoomLock==1
                FixedScale = get(gca, {'xlim', 'ylim'});
            end
            redraw;
        end
    end

    function ratio_change(source, eventdata)
        AX_t=get(axcontrol,'String');
        AX= AX_t{ get(axcontrol,'Value') };
        redraw;
    end

    function cmap_change(source, eventdata)
        CM_t=get(c_map,'String');
        CM= CM_t{ get(c_map,'Value') };
        redraw;
    end

    function redraw
        if scale(1)==scale(2)
            scale(1) = scale(1) - 1;
            scale(2) = scale(2) + 1;
        end
        
        hh=imagesc(I, scale); axis(AX); colormap(CM); colorbar
        
        if zoomLock==1
            set(gca, {'xlim', 'ylim'}, FixedScale);
        end
    end

    function ri_change(source, eventdata)
        ri_t=get(ri_tag,'String');
        ri= ri_t{get(ri_tag,'Value')};
        
        switch ri
            case 'absolute'
                I=abs(Ibuffer);

            case 'phase'
                I=angle(Ibuffer);
                
            case 'real'
                I=real(Ibuffer);
                
            case 'imaginary'
                I=imag(Ibuffer);
        end     
        if zoomLock==1
            FixedScale = get(gca, {'xlim', 'ylim'});
        end
        if scaleLock==0
            
            scaleChange(I);
        end
        redraw;
    end

    function next_pressed(source, eventdata)
        
        if frameN<MaxFrame
            if strcmp(type_text,'matrix')
                frameN=frameN+1;
                Ibuffer=Ioriginal(:,:,frameN);

            else
                frameN=frameN+1;
                Ibuffer=Ioriginal{frameN};
                
            end
            
            if ~isreal(Ibuffer)
                set(ri_tag,'Visible','On');
                ri_change
            else
                set(ri_tag,'Visible','Off');
                I=Ibuffer;
            end
            
            if scaleLock==0
                scaleChange(I);
            end
            if zoomLock==1
                FixedScale = get(gca, {'xlim', 'ylim'});
            end
            redraw;
            set(frame_number,'String',sprintf('Frame #: %d / %d',frameN, MaxFrame));
            set(sld, 'Value', frameN*(1/sno));
        end
        
    end

    function prev_pressed(source, eventdata)
        
        if frameN>1
            if strcmp(type_text,'matrix')
                frameN=frameN-1;
                Ibuffer=Ioriginal(:,:,frameN);

            else
                frameN=frameN-1;
                Ibuffer=Ioriginal{frameN};
            end
            
            if ~isreal(Ibuffer)
                set(ri_tag,'Visible','On');
                ri_change
            else
                set(ri_tag,'Visible','Off');
                I=Ibuffer;
            end
            
            if scaleLock==0
                scaleChange(I);
            end
            if zoomLock==1
                FixedScale = get(gca, {'xlim', 'ylim'});
            end
            redraw;
            set(frame_number,'String',sprintf('Frame #: %d / %d',frameN,MaxFrame));
            set(sld, 'Value', frameN*(1/sno));
        end        
    end

    function scaleChange(img)
        
        if scaleLock==0
            scale = setscale(img);
            slH=uicontrol('style','slider','units','normalized','Position',[0.87 0.1 0.03 0.8],'Max',scale(2), 'Min', scale(1), 'Value', scale(2), 'Callback', {@high_change});
            slL=uicontrol('style','slider','units','normalized','Position',[0.93 0.1 0.03 0.8],'Max',scale(2), 'Min', scale(1), 'Value', scale(1), 'Callback', {@low_change} );
            set(slH,'Max',scale(2),'Min',scale(1),'Value',scale(2));
            set(slL,'Max',scale(2),'Min',scale(1),'Value',scale(1));
        end
        
    end

    function sc_Lock(source, eventdata)
        
        if get(scLock,'Value')
            scaleLock=1;
            set(slH,'Enable','off');
            set(slL,'Enable','off');
        else
            scaleLock=0;
            set(slH,'Enable','on');
            set(slL,'Enable','on');
            scaleChange(I);
            redraw;
        end
    end

    function zoom_Lock(source, eventdata)
        
        if get(zLock,'Value')
            zoomLock=1;
            zoom off  
            set(zoom_handles,'enable','off');
        else
            zoomLock=0;
            reset(gca);
            set(zoom_handles,'enable','on');

            redraw;
        end
    end

    function new_scale=setscale( II ) % Designate current Input image scale as 'scale' parameter.
        
        new_scale(1) = min( II(:) );
        new_scale(2) = max( II(:) );
            
        if new_scale(1) == new_scale(2)
            new_scale(1)= new_scale(1) -1;
            new_scale(2)= new_scale(2) +1;
        end        
        
    end

    function frame_shift(source, event_data)
        
        CurrentFrame=round(get(sld,'value')*MaxFrame);
        if CurrentFrame<=0
            CurrentFrame = 1;
        end
        if strcmp(type_text,'matrix')
            frameN=CurrentFrame;
            Ibuffer=Ioriginal(:,:,frameN);
        else
            frameN=CurrentFrame;
            Ibuffer=Ioriginal{frameN};
        end
        if ~isreal(Ibuffer)
            set(ri_tag,'Visible','On');
            ri_change
        else
            set(ri_tag,'Visible','Off');
            I=Ibuffer;
        end
        if scaleLock==0
            scaleChange(I);
        end
        if zoomLock==1
            FixedScale = get(gca, {'xlim', 'ylim'});
        end
        redraw;
        set(frame_number,'String',sprintf('Frame #: %d / %d',frameN, MaxFrame));

    end

    function showImage(src,evnt)
        if evnt.VerticalScrollCount < 0
            % Scroll Up
            if strcmp(type_text,'matrix')
                frameN=frameN+1;
                if frameN > sno
                	frameN = 1;
                end
                Ibuffer=Ioriginal(:,:,frameN);
            else
                frameN=frameN+1;
                if frameN > sno
                	frameN = 1;
                end
                Ibuffer=Ioriginal{frameN};
            end
            
            if ~isreal(Ibuffer)
                set(ri_tag,'Visible','On');
                ri_change
            else
                set(ri_tag,'Visible','Off');
                I=Ibuffer;
            end
            
        elseif evnt.VerticalScrollCount > 0
            % Scroll Down
            if strcmp(type_text,'matrix')
                frameN=frameN-1;
                if frameN < 1 
                	frameN = sno;
                end
                Ibuffer=Ioriginal(:,:,frameN);

            else
                frameN=frameN-1;
                if frameN < 1 
                	frameN = sno;
                end
                Ibuffer=Ioriginal{frameN};
            end
            
            if ~isreal(Ibuffer)
                set(ri_tag,'Visible','On');
                ri_change
            else
                set(ri_tag,'Visible','Off');
                I=Ibuffer;
            end
        end
        if scaleLock==0
            scaleChange(I);
        end
        if zoomLock==1
            FixedScale = get(gca, {'xlim', 'ylim'});
        end
        redraw;
        set(frame_number,'String',sprintf('Frame #: %d / %d',frameN, MaxFrame) );
        set(sld, 'Value', frameN*(1/(sno)));
    end
 
clim = scale;
cax = ancestor(hh,'axes');

if ~isempty(clim)
  set(cax,'CLim',clim)
elseif ~ishold(cax)
  set(cax,'CLimMode','auto')
end

if nargout > 0
    hhh = hh;
end

end