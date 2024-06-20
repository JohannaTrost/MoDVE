%Create species matrices
clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Parameters that need to be specified/checked before running this script

%Name of the species pool
NameSpeciesPoolSummary='SP_Random_IA_2_IR_60_TimeS_200'; %Name of large species pool, which included population growth rates
NameSpeciesPoolOrig='SP_Random_IntAgeMat_2_IntRec_60_TraitCorrOn'; %Name of original species pool, from which the large species pool was created

%NameSpeciesPool='SpeciesPool_Summary'; %Give meaningful name (the species type is automatically added to the name)
%Name of the specific forest model (under which the single species
%population growth rates are saved in the species pool folder
NameForestModel='ForestModel_Best_50x50';

%The following option defines the type of species pool that is generated:
%0: random species pool (species are randomly selected from the entire trait space)
%1: sequential species pool (species are sequentially selected along a define trait axis)
%2: neutral species pool (all species have the same (defined) traits
SpeciesPoolType=0;

%Define number of species in species pool and total number of species pools to be created
numSpeciesPools=10;
NumberOfSpecies=100;

%Define criteria to select species
NameCriterium='Pop1020-1030';

CriteriaType=1; %1:PopGrowthRate, 2:NumberIndividuals

PopGrowthMin=1.02;
PopGrowthMax=1.03;

IndMeanRelMin=0.9;
IndMeanRelMax=1.1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MainSpeciesPoolFolder='\\OMNISCIENTIA\Members\Gunnar\EpiphyteModels\SpeciesPools';

FileSpeciesPoolGrowthRate=strcat(MainSpeciesPoolFolder,'\',NameSpeciesPoolSummary,'\SpeciesPool_Summary.csv');
SpeciesPool=dlmread(FileSpeciesPoolGrowthRate, '\t');


FileSpeciesPoolHeaders=strcat(MainSpeciesPoolFolder,'\',NameSpeciesPoolSummary,'\SpeciesPool_Headers.xls');
[Test, ColumnHeadersSpeciesPool]=xlsread(FileSpeciesPoolHeaders);

ColPopGrowthMean=find(ismember(ColumnHeadersSpeciesPool,'PopGrowth_Mean'));
ColIndMeanRel=find(ismember(ColumnHeadersSpeciesPool,'IndMeanRel_Mean'));

%Make subset containing only species fullfilling the criteria
if(CriteriaType==1)
    SpeciesPoolSelected=SpeciesPool(SpeciesPool(:,ColPopGrowthMean)>=PopGrowthMin & SpeciesPool(:,ColPopGrowthMean)<=PopGrowthMax,:);
elseif(CriteriaType==2)
    SpeciesPoolSelected=SpeciesPool(SpeciesPool(:,ColIndMeanRel)>=IndMeanRelMin & SpeciesPool(:,ColIndMeanRel)<=IndMeanRelMin,:);
end

%Test if the number of species fullfilling the criteria larger than the
%desired number of species
NumberSpeciesCriteria=size(SpeciesPoolSelected,1);
if(NumberSpeciesCriteria<NumberOfSpecies)
    disp('Not enough species fullfilling the criteria => rework!!')
    return
end

%Create general folder for specific species pools if it does not exist
SaveDirectoryGeneral=strcat(MainSpeciesPoolFolder,'\',NameSpeciesPoolSummary);
if exist(SaveDirectoryGeneral, 'dir')==0
    mkdir(SaveDirectoryGeneral);
end

SaveDirectory=strcat(SaveDirectoryGeneral,'\',NameForestModel,'_',NameCriterium);
if exist(SaveDirectory, 'dir')==0
    mkdir(SaveDirectory);
end


%Maximum colum of regular species pool
ColMax=find(ismember(ColumnHeadersSpeciesPool,'AgeAtMaturity'));

%Select species randomly and save in specific species pool folder
for Num=1:numSpeciesPools

    %Randomly select species from species pool
    RandomSpecies=randsample(NumberSpeciesCriteria,NumberOfSpecies);
    SpeciesPoolRandDetailed=SpeciesPoolSelected(RandomSpecies,:);
    SpeciesPoolRand=SpeciesPoolSelected(RandomSpecies,1:ColMax);
    SpeciesPoolRandDetailed(:,1)=1:NumberOfSpecies;
    SpeciesPoolRand(:,1)=1:NumberOfSpecies;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Save species pools/trait matrices
    %Create dataset from table (including headers)
    SpeciesTable=dataset({SpeciesPoolRand,ColumnHeadersSpeciesPool{1:ColMax}});
    SpeciesTableDetailed=dataset({SpeciesPoolRandDetailed,ColumnHeadersSpeciesPool});

    %Save species pools/trait matrix
    SaveFile=strcat(SaveDirectory,'\SpeciesPool',num2str(Num),'.csv');
    dlmwrite(SaveFile,SpeciesPoolRand,'\t');

    SaveFileDataset=strcat(SaveDirectory,'\SpeciesPoolHeader',num2str(Num),'.txt');
    export(SpeciesTable,'file',SaveFileDataset)

    SaveFileDatasetDetailed=strcat(SaveDirectory,'\SpeciesPoolDetailedHeader',num2str(Num),'.txt');
    export(SpeciesTableDetailed,'file',SaveFileDatasetDetailed)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

%Save ColumnHeader for further use
xlswrite(strcat(SaveDirectory,'\ColumnHeaders.xls'),ColumnHeadersSpeciesPool(1:ColMax));
xlswrite(strcat(SaveDirectory,'\ColumnHeadersDetailed.xls'),ColumnHeadersSpeciesPool);

%Save SpeciesPoolType for further use
save(strcat(SaveDirectory,'\SpeciesPoolType.mat'),'SpeciesPoolType')

%Save number of species pools and number of species for further use
save(strcat(SaveDirectory,'\numSpeciesPools.mat'),'numSpeciesPools')
save(strcat(SaveDirectory,'\NumberOfSpecies.mat'),'NumberOfSpecies')

copyfile(strcat(MainSpeciesPoolFolder,'\',NameSpeciesPoolOrig,'\TraitRanges.csv'),...
    strcat(SaveDirectory,'\TraitRanges.csv'));
