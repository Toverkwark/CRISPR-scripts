function ProcessScreen(ReadCountFile,ConstructStatisticsFile,AggregateStatisticsFile,ControlSamples,TreatedSamples,SampleNames,MCPlotName)
EssentialGenes={'RPS9';'RPS8';'RPS7';'RPS3A';'RPS27';'RPS24';'RPS19';'RPS17';'RPS13';'RPS11';'RPLP1';'RPL9';'RPL6';'RPL5';'RPL36';'RPL35A';'RPL34';'RPL30';'RPL3';'RPL27';'RPL19';'RPL18A';'RPL11';'PSMD7';'PSMD6';'PSMD11';'PSMD1';'PSMC4';'PSMC2';'PSMC1';'PSMB3';'PSMB2';'PSMA3';'POLR2F';'POLR2D';'POLR2A';'POLA1';'NUP98';'NUP93';'NUP54';'NUP205';'NUP133';'KPNB1';'COPZ1';'COPS8';'COPS6';'COPS4';'COPS2';'COPB1';'COPA';};
NonEssentialGenes={'CRYGB';'KRT77';'DMRTB1';'POTEA';'NLRP5';'VN1R5';'OR9Q2';'TAAR8';'OR12D2';'LUZP4';'TGM6';'SAGE1';'TPH2';'LHX5';'TAS2R13';'VN1R2';'DEFB129';'RXFP2';'ADH7';'DMRTC2';'RNASE9';'ABCG8';'PLA2G2E';'KRT74';'IL22';'DPCR1';'TAAR1';'TAS2R9';'CYP7A1';'MAGEB3';'NPSR1';'OLIG2';'MRGPRD';'CABP5';'POU4F2';'OR52E8';'TRIM42';'OC90';'HTR3D';'RPTN';'IL1F10';'LYZL6';'OTUD6A';'KRT25';'KRT9';'FCRL4';'SPATA16';'NPHS2';'FAM71B';'PIWIL3';};
NumberOfSamples=size(ControlSamples,2)+size(TreatedSamples,2);

%Read in sequencing data and sort by library ID
disp('Reading in sequencing data');
fID=fopen(ReadCountFile);
tmp=textscan(fID,'%s%s%f%f%f%f%f%f%f%f%f%f%f%f','HeaderLines',1,'Delimiter','\t');
fclose(fID);
LibraryIDs=tmp{1};
Genes=tmp{2};
ReadCounts=cell2mat(tmp(3:14));
[tmp idx]=sort(LibraryIDs);
LibraryIDs=LibraryIDs(idx);
ReadCounts=ReadCounts(idx,:);
Genes=Genes(idx);

%Select samples
ReadCounts=ReadCounts(:,[ControlSamples TreatedSamples]);
% SampleNames=SampleNames([ControlSamples TreatedSamples]);

%Read in Construct statistics data and sort by library ID
fID=fopen(ConstructStatisticsFile);
tmp=textscan(fID,'%s%s%s%s%f%f%f%f%f%f%f%f%f%s','HeaderLines',1);
LibraryIDsConstructStatistics=tmp{1};
PValuesConstructStatistics=tmp{10};
[tmp idx]=sort(LibraryIDsConstructStatistics);
LibraryIDsConstructStatistics=LibraryIDsConstructStatistics(idx);
PValuesConstructStatistics=PValuesConstructStatistics(idx);
fclose(fID);

%Read in gene statistics data
fID=fopen(AggregateStatisticsFile);
tmp=textscan(fID,'%s%f%f%f%f%f%f%f%f%f%f%f','HeaderLines',1);
NameGeneStatistics=tmp{1};
PValuesGeneStatistics=tmp{4};
fclose(fID);

%Assert that LibraryIDs and LibraryIDsStatistics are identical
disp('Asserting library equality');
for i=1:size(LibraryIDs,1)
    if(~strcmp(LibraryIDs{i},LibraryIDsConstructStatistics{i}) || (size(LibraryIDs,1) ~= size(LibraryIDsConstructStatistics,1)))
        disp('Equal library assertion failed');
        return;
    end
end

%Sort all data according to Construct P value
[tmp, idx]=sort(PValuesConstructStatistics);
LibraryIDs=LibraryIDs(idx);
ReadCounts=ReadCounts(idx,:);
Genes=Genes(idx,:);
PValuesConstructStatistics=PValuesConstructStatistics(idx,:);

%Perform Benjamini Hochberg correction for construct P-values
PValuesConstructStatistics=mafdr(PValuesConstructStatistics,'BHFDR',true);

%Map essentialities
disp('Mapping essentialities');
for i=1:size(LibraryIDs,1)
    if (find(strcmp(Genes{i},EssentialGenes)))
        Essentialities(i)=1;
    else
        if (strfind(LibraryIDs{i},'NonTargetingControl'))
            Essentialities(i)=2;
        else
            Essentialities(i)=0;
        end
    end
end

%Plot correlation graphs
disp('Plotting correlation graphs');
figure('Position',get(0,'ScreenSize'));
for row=1:NumberOfSamples
    for column=1:NumberOfSamples
        subaxis(NumberOfSamples,NumberOfSamples,column+(row-1)*NumberOfSamples,'Spacing',0.025,'Padding',0,'Margin',0.06);
        scatter(log2(ReadCounts(:,row)),log2(ReadCounts(:,column)),9,'.');
        set(gca,'FontSize',6);
        xlabel(['^2Log read counts ' SampleNames{row}]);
        ylabel(['^2Log read counts ' SampleNames{column}]);
        clear CorrData;
        CorrData(:,1)=(ReadCounts(:,row));
        CorrData(:,2)=(ReadCounts(:,column));
        CorrData=CorrData(find(CorrData(:,1)>1),:);
        CorrData=CorrData(find(CorrData(:,2)>1),:);
        CorrelationCoefficient=corrcoef(log2(CorrData(:,1)),log2(CorrData(:,2)));
        text(0.1*log2(max(ReadCounts(:))),0.9*log2(max(ReadCounts(:))),['\rho=' num2str(CorrelationCoefficient(1,2))]);
        set(gca,'XLim',[0 log2(max(ReadCounts(:)))]);
        set(gca,'YLim',[0 log2(max(ReadCounts(:)))]);
    end
end
% close;

%Normalize reads
disp('Normalizing read counts');
NormalizedReads=NormalizeReads(ReadCounts);

%Average replicates
disp('Averaging replicates');
AveragedData(:,1)=geomean(NormalizedReads(:,1:size(ControlSamples,2))')';
AveragedData(:,2)=geomean(NormalizedReads(:,size(ControlSamples,2)+1:size(ControlSamples,2)+size(TreatedSamples,2))')';

%Produce MC plots mapped for essentiality or P value order
disp('Producing MC Plot');
MCPlotSVG(AveragedData(:,2),AveragedData(:,1),LibraryIDs,Genes,Essentialities,[MCPlotName '.html']);

%Sort Construct data based on Construct P Values
disp('Plotting Construct waterfall');
[tmp, idx]=sort(PValuesConstructStatistics);
LibraryIDs=LibraryIDs(idx);
Genes=Genes(idx,:);
Essentialities=Essentialities(idx);
PValuesConstructStatistics=PValuesConstructStatistics(idx);
AveragedData=AveragedData(idx,:);
NormalizedReads=NormalizedReads(idx,:);

%Plot Construct P-Value beeswarm plot
figure;
[GroupIDs GroupNames]=grp2idx(Genes);
for i=1:size(GroupIDs,1)
    GroupIDs(i)=find(strcmp(GroupNames(GroupIDs(i)),NameGeneStatistics));
end
plotSpread(log10(PValuesConstructStatistics),'distributionIdx',GroupIDs,'categoryIdx',Essentialities','categoryColors',[[0 0.6 0.5];[0.9 0.6 0];[0.35 0.7 0.9]]);
XTickLabels={};
for i=50:50:size(PValuesConstructStatistics,1)
    XTickLabels=[XTickLabels;num2str(i)];
end
set(gca,'XTick',50:50:size(PValuesConstructStatistics,1),'XTickLabel',XTickLabels,'FontName','Arial');
xlabel('Gene ID','FontName','Arial');
ylabel('^1^0Log Construct FDR','FontName','Arial');
% Children=get(gca,'Children');
% for i=1:size(Children,1)
%     set(Children(i),'MarkerSize',20)
% end


%Gather ROC info for Constructs
PositivesFound=0;
NegativesFound=0;
ConstructTPR(size(Essentialities,2))=0;
ConstructFPR(size(Essentialities,2))=0;
TotalPositives=sum(Essentialities==1);
TotalNegatives=sum(Essentialities==2);
for i=1:size(Essentialities,2)
    if(Essentialities(i)==1)
        PositivesFound=PositivesFound+1;
    end
    if(Essentialities(i)==2)
        NegativesFound=NegativesFound+1;
    end
    ConstructTPR(i)=PositivesFound/TotalPositives;
    ConstructFPR(i)=NegativesFound/TotalNegatives;
end

%Add 0,0 as a datapoint
ConstructTPR=[0 ConstructTPR];
ConstructFPR=[0 ConstructFPR];
figure;
plot(ConstructFPR,ConstructTPR,'Color',[0.9 0.6 0],'LineWidth',2.9);
xlabel('Construct based FPR');
ylabel('Construct based TPR');
title('Construct based ROC curve');

%Plot gene waterfall plot colored by essentiality
disp('Plotting gene waterfall');
for i=1:size(NameGeneStatistics,1)
    if(find(strcmp(NameGeneStatistics{i},EssentialGenes)))
        GeneEssentialityStatistics(i)=1;
    else
        if(find(strcmp(NameGeneStatistics{i},NonEssentialGenes)))
            GeneEssentialityStatistics(i)=2;
        else
            if(strcmp(NameGeneStatistics{i},'NonTargetingControl'))
                GeneEssentialityStatistics(i)=2;
            else
                GeneEssentialityStatistics(i)=0;
            end
        end
    end
end

figure;
%Perform Benjamini Hochberg correction for gene P-values
PValuesGeneStatistics=mafdr(PValuesGeneStatistics,'BHFDR',true);
EssentialAggregates=log10(PValuesGeneStatistics);
NonEssentialAggregates=log10(PValuesGeneStatistics);
OtherAggregates=log10(PValuesGeneStatistics);
for i=1:size(PValuesGeneStatistics,1)
    if(GeneEssentialityStatistics(i)==1)
        NonEssentialAggregates(i)=0;
        OtherAggregates(i)=0;
    else
        if(GeneEssentialityStatistics(i)==2)
            EssentialAggregates(i)=0;
            OtherAggregates(i)=0;
        else
            EssentialAggregates(i)=0;
            NonEssentialAggregates(i)=0;
        end
    end
end
h=bar(EssentialAggregates);
set(h,'FaceColor',[0.9 0.6 0]);
set(h,'LineStyle','none');
set(h,'BarWidth',1);
hold on;
h=bar(NonEssentialAggregates);
set(h,'FaceColor',[0.35 0.7 0.9]);
set(h,'LineStyle','none');
set(h,'BarWidth',1);
h=bar(OtherAggregates);
set(h,'FaceColor',[0 0.6 0.5]);
set(h,'LineStyle','none');
set(h,'BarWidth',1);
ylabel('^1^0Log Gene FDR','FontName','Arial');
xlabel('Gene ID','FontName','Arial');
set(gca,'YLim',[-5 0]);
hold on;
scatter(1:size(GeneEssentialityStatistics,2),zeros(size(GeneEssentialityStatistics,2),1),50,GeneEssentialityStatistics,'.')
colormap([[0 0.6 0.5];[0.9 0.6 0];[0.35 0.7 0.9]]);
set(gca,'FontName','Arial');

%Plot beeswarm plots of fold change distributions 
for i=1:size(NameGeneStatistics,1)
    BeeswarmData{:,i}=log2(AveragedData(find(strcmp(NameGeneStatistics{i},Genes)),2)./AveragedData(find(strcmp(NameGeneStatistics{i},Genes)),1));
end
figure;
plotSpread(BeeswarmData,'distributionIdx',GroupIDs,'categoryIdx',Essentialities','categoryColors',[[0 0.6 0.5];[0.9 0.6 0];[0.35 0.7 0.9]]);
set(gca,'XTick',50:50:size(NameGeneStatistics,1));
set(gca,'YLim',[-12 4]);
line([0 size(NameGeneStatistics,1)],[0 0],'Color',[0 0 0]);
ylabel('^2Log fold change');
xlabel('Gene ID');
