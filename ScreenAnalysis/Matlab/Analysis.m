clear all;
close all;
ReadCountFile='..\SequencingData\crispr.epigenome.output';
ConstructStatisticsFile='..\MageckAnalysis\t1_vs_t0_2D.sgrna_summary.txt';
AggregateStatisticsFile='..\MageckAnalysis\t1_vs_t0_2D.gene_summary.txt';
ControlSamples=[1 2 3];
TreatedSamples=[4 5 6];
SampleNames={   't=0 Rep.1',
                't=0 Rep.2',
                't=0 Rep.3',
                't=1 Rep.1 2D',
                't=1 Rep.2 2D',
                't=1 Rep.3 2D'
};
ScreenTitle='t1_vs_t0_2D';
ProcessScreen(ReadCountFile,ConstructStatisticsFile,AggregateStatisticsFile,ControlSamples,TreatedSamples,SampleNames,ScreenTitle);
