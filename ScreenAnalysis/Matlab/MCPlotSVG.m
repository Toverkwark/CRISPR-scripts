function MCPlotSVG(Treated,Control,LibraryIDs,Genes,Colors,FileName)
HorizontalPixels=666;
VerticalPixels=500;
MinimalX=4;
MaximalX=12;
if(findstr('UMUC',FileName))
    MaximalX=16;
end
MinimalY=-3;
MaximalY=3;
PixelSize=2;
XOffset=10;
YInterval=0.2;
HistogramOffset=10;
HistogramWidth=150;

%Add 1 to everything to prevent log problems
Treated=Treated+1;
Control=Control+1;
A=log2(Control);
M=-log2(Treated./Control);
LineFeed=char(10);

%Create header
clear OutputText;
fID=fopen('SVGHEADER.txt');
Line=fgetl(fID);
i=0;
while ischar(Line)
    i=i+1;
    Header{i}=[Line LineFeed];
    Line=fgetl(fID);
end
fclose(fID);
Header=Header';
%Create axes
Header=[Header;'<line x1="' num2str(XOffset) '" y1="0" x2="' num2str(HorizontalPixels+XOffset)  '" y2="0" stroke-width="0.5" stroke="black" />' LineFeed];
Header=[Header;'<line x1="' num2str(XOffset) '" y1="0" x2="' num2str(XOffset) '" y2="' num2str(VerticalPixels) '" stroke-width="0.5" stroke="black" />' LineFeed];
Header=[Header;'<line x1="' num2str(HorizontalPixels+XOffset) '" y1="0" x2="' num2str(HorizontalPixels+XOffset) '" y2="' num2str(VerticalPixels) '" stroke-width="0.5" stroke="black" />' LineFeed];
Header=[Header;'<line x1="' num2str(XOffset) '" y1="' num2str(VerticalPixels) '" x2="' num2str(HorizontalPixels+XOffset) '" y2="' num2str(VerticalPixels) '" stroke-width="0.5" stroke="black" />' LineFeed];
for i=MinimalY:MaximalY 
	VerticalLevel=(-i-MinimalY)*VerticalPixels/(MaximalY-MinimalY);
	Header=[Header;'<line x1="' num2str(XOffset) '" y1="' num2str(VerticalLevel) '" x2="' num2str(HorizontalPixels+XOffset) '" y2="' num2str(VerticalLevel) '" stroke-width=".2" stroke="black" />' LineFeed];
	YAxisLabel=num2str(i);
    if(i==MinimalY)
        YAxisLabel=['&#8804;' num2str(i)];
    else
        if(i==MaximalY)
            YAxisLabel=['&#8805;' num2str(i)];
        else
            YAxisLabel=num2str(i);
        end
    end
    Header=[Header;'<text font-family="Arial" font-size="8" x="5" y="' num2str(VerticalLevel+0.005*VerticalPixels) '" text-anchor="end">' YAxisLabel '</text>' LineFeed];
    if i==0
        Header=[Header;'<line x1="' num2str(XOffset) '" y1="' num2str(VerticalLevel) '" x2="' num2str(HorizontalPixels+XOffset) '" y2="' num2str(VerticalLevel) '" stroke-width="2" stroke="black" />' LineFeed];
    end
end
for i=MinimalX:MaximalX
	HorizontalLevel=XOffset+((i-MinimalX)*HorizontalPixels/(MaximalX-MinimalX));
	Header=[Header;'<line y1="0" x1="' num2str(HorizontalLevel) '" y2="' num2str(VerticalPixels) '" x2="' num2str(HorizontalLevel) '" stroke-width=".2" stroke="black" />' LineFeed];
    XAxisLabel=num2str(i);
    if(i==MinimalX)
        XAxisLabel=['&#8804;' num2str(i)];
    else
        if(i==MaximalX)
            XAxisLabel=['&#8805;' num2str(i)];
        else
            XAxisLabel=num2str(i);
        end
    end
	Header=[Header;'<text font-family="Arial" font-size="8" x="' num2str(HorizontalLevel) '" y="' num2str(-0.005*VerticalPixels) '" text-anchor="middle">' XAxisLabel '</text>' LineFeed];
end
%Write dots
XPosition=(A-MinimalX).*(HorizontalPixels/(MaximalX-MinimalX));
YPosition=(M-MinimalY).*(VerticalPixels/(MaximalY-MinimalY));
XPosition=max(XPosition,0);
XPosition=min(XPosition,HorizontalPixels);
YPosition=max(YPosition,0);
YPosition=min(YPosition,VerticalPixels);

%Map colors if necessary
CurrentColormap=[[0 0.6 0.5];[0.9 0.6 0];[0.35 0.7 0.9]];
Colors=Colors';
% CurrentColormap=jet;
if(size(Colors,2)==1)
    MinColor=min(Colors);
    MaxColor=max(Colors);
    for i=1:size(Colors,1)
        CurrentColor=Colors(i);
        CurrentColor=(CurrentColor-MinColor)/MaxColor;
        CurrentColor=CurrentColor*(size(CurrentColormap,1)-1);
        CurrentColor=1+floor(CurrentColor+0.5);
        MappedColors(i,:)=CurrentColormap(CurrentColor,:);
    end
else
    MappedColors=Colors;
end

%Write all the dots
OutputText{size(A,1)}='';
for i=1:size(A,1)
    FillColor=floor(MappedColors(i,:).*255+0.5);
    %Comment out this if end statement for only plotting controls
%     if(sum(MappedColors(i,:)~=CurrentColormap(1,:))>0)
        OutputText{i}=['<circle class="' Genes{i} '" pixelsize="' num2str(PixelSize) '" libraryid="' LibraryIDs{i} '" selected="no" cx="' num2str(XPosition(i)+XOffset) '" cy="' num2str(YPosition(i)) '" r="' num2str(PixelSize) '" stroke="black" stroke-width=".1" fill="rgb(' num2str(FillColor(1)) ',' num2str(FillColor(2)) ',' num2str(FillColor(3)) ')" onclick="ClickDetected(evt)" onmouseover="MouseOverDetected(evt)" onmouseout="MouseOutDetected(evt)"/>' LineFeed];
%     end
end

%Plot histogram
[n1,xout]=hist(M(find(Colors==0)),MinimalY:YInterval:MaximalY);
[n2,xout]=hist(M(find(Colors==1)),MinimalY:YInterval:MaximalY);
[n3,xout]=hist(M(find(Colors==2)),MinimalY:YInterval:MaximalY);
% binwidths = diff([MinimalY xout(1:end-1)+diff(xout)/2 MaximalY]);
n1 = n1./sum (n1);
n2 = n2./sum (n2);
n3 = n3./sum (n3);

% MaxPeak=max(max(n1),max(n2));
Path1='<path d="';
for i=1:size(xout,2)
    PlotYPosition=(xout(i)-MinimalY).*(VerticalPixels/(MaximalY-MinimalY));
    PlotXPosition=XOffset+HorizontalPixels+HistogramOffset+(HistogramWidth*n1(i));
    FillColor=floor(CurrentColormap(1,:).*255+0.5);
    if (i==1) 
        Path1=[Path1 'M ' num2str(PlotXPosition) ' ' num2str(PlotYPosition) ' '];
    else
        Path1=[Path1 'L ' num2str(PlotXPosition) ' ' num2str(PlotYPosition) ' '];
    end
end
Path1=[Path1 '" stroke="rgb(' num2str(FillColor(1)) ',' num2str(FillColor(2)) ',' num2str(FillColor(3)) ')" stroke-width="2" fill="none"/>' LineFeed];
Path2='<path d="';
for i=1:size(xout,2)
    PlotYPosition=(xout(i)-MinimalY).*(VerticalPixels/(MaximalY-MinimalY));
    PlotXPosition=XOffset+HorizontalPixels+HistogramOffset+(HistogramWidth*n2(i));
    FillColor=floor(CurrentColormap(2,:).*255+0.5);
    if (i==1) 
        Path2=[Path2 'M ' num2str(PlotXPosition) ' ' num2str(PlotYPosition) ' '];
    else
        Path2=[Path2 'L ' num2str(PlotXPosition) ' ' num2str(PlotYPosition) ' '];
    end
end
Path2=[Path2 '" stroke="rgb(' num2str(FillColor(1)) ',' num2str(FillColor(2)) ',' num2str(FillColor(3)) ')" stroke-width="2" fill="none"/>' LineFeed];
Path3='<path d="';
for i=1:size(xout,2)
    PlotYPosition=(xout(i)-MinimalY).*(VerticalPixels/(MaximalY-MinimalY));
    PlotXPosition=XOffset+HorizontalPixels+HistogramOffset+(HistogramWidth*n3(i));
    FillColor=floor(CurrentColormap(3,:).*255+0.5);
    if (i==1) 
        Path3=[Path3 'M ' num2str(PlotXPosition) ' ' num2str(PlotYPosition) ' '];
    else
        Path3=[Path3 'L ' num2str(PlotXPosition) ' ' num2str(PlotYPosition) ' '];
    end
end
Path3=[Path3 '" stroke="rgb(' num2str(FillColor(1)) ',' num2str(FillColor(2)) ',' num2str(FillColor(3)) ')" stroke-width="2" fill="none"/>' LineFeed];

OutputText{size(OutputText,2)+1}=Path3;
OutputText{size(OutputText,2)+1}=Path2;
OutputText{size(OutputText,2)+1}=Path1;

%Write outputfile
OutputText=OutputText';
fID=fopen(FileName,'w+');
for i=1:size(Header,1)
    fwrite(fID,Header{i});
end
for i=1:size(OutputText,1)
    fwrite(fID,OutputText{i});
end
fwrite(fID,['</svg>' LineFeed]);
fclose(fID);

end
