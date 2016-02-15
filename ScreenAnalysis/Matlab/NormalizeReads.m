function NormalizedReads=NormalizeReads(RawReads)
%First, calculate a reference sample by taking the geomean of all non-zero
%values
for i=1:size(RawReads,1)
    Reference(i)=geomean(RawReads(i,RawReads(i,:)>0));
    if(isnan(Reference(i)))
        Reference(i)=1;
    end
end

%Prevent long processing times by pre allocating vectors
SizeFactors(size(RawReads,1),size(RawReads,2))=0;
NormalizedData(size(RawReads,1),size(RawReads,2))=0;

%Scale all data to the reference sample
for i=1:size(RawReads,1)
    SizeFactors(i,:)=RawReads(i,:)./Reference(i);
end

%Determine the SizeFactors, which is the median of all non-zero data 
for i=1:size(SizeFactors,2)
    SizeFactor(i)=median(SizeFactors(SizeFactors(:,i)>0,i));
end
SizeFactor=1./SizeFactor';

%Normalize the data
for i=1:size(RawReads,2)
    NormalizedReads(:,i)=RawReads(:,i).*SizeFactor(i);
end

