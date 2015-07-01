function MCPlotSVG(Treated,Control,LibraryIDs,Genes,Colors,FileName)
HorizontalPixels=666;
VerticalPixels=500;
MinimalX=0;
MaximalX=15;
MinimalY=-10;
MaximalY=10;
PixelSize=2;

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
    OutputText{i}=[Line LineFeed];
    Line=fgetl(fID);
end
fclose(fID);
OutputText=OutputText';

%Create axes
OutputText=[OutputText;'<line x1="0" y1="0" x2="' num2str(HorizontalPixels)  '" y2="0" stroke-width="0.5" stroke="black" />' LineFeed];
OutputText=[OutputText;'<line x1="0" y1="0" x2="0" y2="' num2str(VerticalPixels) '" stroke-width="0.5" stroke="black" />' LineFeed];
OutputText=[OutputText;'<line x1="' num2str(HorizontalPixels) '" y1="0" x2="' num2str(HorizontalPixels) '" y2="' num2str(VerticalPixels) '" stroke-width="0.5" stroke="black" />' LineFeed];
OutputText=[OutputText;'<line x1="0" y1="' num2str(VerticalPixels) '" x2="' num2str(HorizontalPixels) '" y2="' num2str(VerticalPixels) '" stroke-width="0.5" stroke="black" />' LineFeed];
for i=MinimalY:MaximalY-1 
	VerticalLevel=(-i-MinimalY)*VerticalPixels/(MaximalY-MinimalY);
	OutputText=[OutputText;'<line x1="0" y1="' num2str(VerticalLevel) '" x2="' num2str(HorizontalPixels) '" y2="' num2str(VerticalLevel) '" stroke-width=".2" stroke="black" />' LineFeed];
	OutputText=[OutputText;'<text font-family="Arial" font-size="8" x="' num2str(0.005*HorizontalPixels) '" y="' num2str(VerticalLevel-0.005*VerticalPixels) '">' num2str(i) '</text>' LineFeed];
    if i==0
        OutputText=[OutputText;'<line x1="0" y1="' num2str(VerticalLevel) '" x2="' num2str(HorizontalPixels) '" y2="' num2str(VerticalLevel) '" stroke-width="2" stroke="black" />' LineFeed];
    end
end
for i=MinimalX+1:MaximalX
	HorizontalLevel=(i-MinimalX)*HorizontalPixels/(MaximalX-MinimalX);
	OutputText=[OutputText;'<line y1="0" x1="' num2str(HorizontalLevel) '" y2="' num2str(VerticalPixels) '" x2="' num2str(HorizontalLevel) '" stroke-width=".2" stroke="black" />' LineFeed];
	OutputText=[OutputText;'<text font-family="Arial" font-size="8" x="' num2str(HorizontalLevel-0.005*HorizontalPixels) '" y="' num2str(-0.005*VerticalPixels) '">' num2str(i) '</text>' LineFeed];
end
%Write dots
XPosition=(A-MinimalX).*(HorizontalPixels/(MaximalX-MinimalX));
YPosition=(M-MinimalY).*(VerticalPixels/(MaximalY-MinimalY));
XPosition=max(XPosition,0);
XPosition=min(XPosition,HorizontalPixels);
YPosition=max(YPosition,0);
YPosition=min(YPosition,VerticalPixels);

%Map colors if necessary
CurrentColormap=[[0 0 1]; [1 0 0]; [0 1 0]];
% Colors=Colors';
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

for i=1:size(A,1)
    FillColor=floor(MappedColors(i,:).*255+0.5);
	OutputText=[OutputText;'<circle class="' Genes{i} '" pixelsize="' num2str(PixelSize) '" libraryid="' LibraryIDs{i} '" selected="no" cx="' num2str(XPosition(i)) '" cy="' num2str(YPosition(i)) '" r="' num2str(PixelSize) '" stroke="black" stroke-width=".1" fill="rgb(' num2str(FillColor(1)) ',' num2str(FillColor(2)) ',' num2str(FillColor(3)) ')" onclick="ClickDetected(evt)" onmouseover="MouseOverDetected(evt)" onmouseout="MouseOutDetected(evt)"/>' LineFeed];
end

%Write outputfile
fID=fopen(FileName,'w+');
for i=1:size(OutputText,1)
    fwrite(fID,OutputText{i});
end
fwrite(fID,['</svg>' LineFeed]);
fclose(fID);

end
