function shape=readoff(fname,displayon)
% r kakarala
if (nargin < 1 || isempty(fname))
    % rad a shape in .off file
    [fname,pname]=uigetfile('*.off','Get shape in .off file');
    fname = [pname fname];
end;
if nargin < 2
    displayon = 1;
end;

%open the off file
dfile=fopen(fname,'r');

%file starter 'OFF'
Lstr=fgetl(dfile); %#ok<*NASGU>

%3D model info
Lstr=fgetl(dfile);
Ldata=str2num(Lstr); %#ok<*ST2NM>
shape.info=Ldata;

%readin vertice
for j=1:shape.info(1)
    Lstr=fgetl(dfile);
    Ldata=str2num(Lstr);
    shape.ver(j,:)=Ldata;
end

%readin faces
for j=1:shape.info(2)
    Lstr=fgetl(dfile);
    Ldata=str2num(Lstr);
    shape.fac(j,:)=Ldata;
end

fclose(dfile);

if (displayon)
    
    handle = patch ( 'Vertices', shape.ver, 'Faces', shape.fac+1 );
    % shape.fac + 1 because off files indexes faces starting from 0
    set ( handle, 'FaceColor', [0.5, 0.6, 0.8], 'EdgeColor', 'Black' );
    
    axis equal;
    axis off;
    view(3);
    %grid on;
    
    xlabel ( '--X axis--' )
    ylabel ( '--Y axis--' )
    zlabel ( '--Z axis--' )
    
    title(fname);
    
end;