%Create the initial epiphyte distrubution depending on the epiphyte traits and the initial microhabitat matrix
clear all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Parameters that need to be specified/checked before running this script

%Name of epiphyte model 
%FolderEpiphyteModel='EM_20160213'; %I should think about naming
SingleSpeciesModel=0; %1: Single species model, 0: Community model

%Name of species pool 
FolderSpeciesPools='SP_Random_IntAgeMat_2_IntRec_60_TraitCorrOn';

%Directory where model is save and directory where microhabitat matrices
%are stored
DirectoryModelMain='C:\Gunnar\EpiphyteModels\EpiphyteModels';
DirectoryMicrohabitatMain='C:\Gunnar\EpiphyteModels\MicrohabitatMatrices';
DirectorySpeciesPoolsMain='\\OMNISCIENTIA\Members\Gunnar\EpiphyteModels\SpeciesPools';

%Name of microhabitat matrix 
%FolderMicrohabitat='Microhabitat_ForestModel_20160113_50x50x80_timesteps_200-1000';
FolderMicrohabitat='ForestModel_Best_50x50';
Replicate=4;
MicrohabitatType=1; %Define which type of forest the microhabitat belongs to. 1: dynamic forest, 2: static forest, 3: uniform forest

%Choose species pools to use and number of replicates per species pool
numSpeciesPools=[91,100]; %Start and end number of  species pools
replicatePerSpeciesPool=5; %Number of replicates per species pool
TimeStep=200; %Time step for which the Initial distribution is generated

%The suitable voxel can either be the voxel
%with the highest available surface area (MethodVoxel=1), or a random voxel
%MethodVoxel=0)
MethodVoxel=0;

%Define how many individuals per species are used, and how many of them are initially mature
%This variable defines if the NumberSpecies are total numbers irrespective
%of the model area (ScalingPerHa=0), or if the NumberSpecies or given per
%hectar and are scaled to the model area (ScalingPerHa=1)
ScalingPerHa=0; 
IndividualsPerSpecies=100;
PercentageMaturePerSpecies=50;

%This parameter set the scaling between the 
SurfaceBiomassScaling=10000*0.01; %cm^2 per m^2
Imax=900;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Folders and directories (these files should not change)

%The maximum path name length can be exeeded thus the folder names are abbreviated
if MicrohabitatType==1
   FolderEpiphyteModel='DynamicForests';
elseif MicrohabitatType==2
   FolderEpiphyteModel='StaticForests';
elseif MicrohabitatType==3
   FolderEpiphyteModel='UniformForests';
end

%Create main model folder under which the initial distribution is saved
if SingleSpeciesModel==1
    DirectoryEpiphyteModelMain1=strcat(DirectoryModelMain,'\SingleSpeciesModels');
    DirectoryEpiphyteModelMain2=strcat(DirectoryModelMain,'\SingleSpeciesModels\',FolderEpiphyteModel);
    DirectoryEpiphyteModel=strcat(DirectoryModelMain,'\SingleSpeciesModels\',FolderEpiphyteModel,'\',FolderMicrohabitat,'_Rep',num2str(Replicate));
elseif SingleSpeciesModel==0
    DirectoryEpiphyteModelMain1=strcat(DirectoryModelMain,'\CommunityModels');
    DirectoryEpiphyteModelMain2=strcat(DirectoryModelMain,'\CommunityModels\',FolderEpiphyteModel);
    DirectoryEpiphyteModel=strcat(DirectoryModelMain,'\CommunityModels\',FolderEpiphyteModel,'\',FolderMicrohabitat,'_Rep',num2str(Replicate));
end

%Generate abbreveation for species pool (in IniDist folder)
PosUnderscores=regexp(FolderSpeciesPools, '_');
NameSpeciesPoolSave=strcat('SP_',FolderSpeciesPools((PosUnderscores(1)+1):(PosUnderscores(2)-1)),...
    '_IA_',FolderSpeciesPools((PosUnderscores(3)+1):(PosUnderscores(4)-1)),...
    '_IR_',FolderSpeciesPools((PosUnderscores(5)+1):(PosUnderscores(6)-1)),...
    '_TimeS_',num2str(TimeStep));

DirectoryIntitalDistributionMain=strcat(DirectoryEpiphyteModel,'\IniDist');
DirectoryIntitalDistribution=strcat(DirectoryIntitalDistributionMain,'\',NameSpeciesPoolSave);

%Generate directories to save the model if not existent
if exist(DirectoryEpiphyteModelMain1, 'file')==0
    mkdir(DirectoryEpiphyteModelMain1); end
if exist(DirectoryEpiphyteModelMain2, 'file')==0
    mkdir(DirectoryEpiphyteModelMain2); end
if exist(DirectoryEpiphyteModel, 'file')==0
    mkdir(DirectoryEpiphyteModel); end
if exist(DirectoryIntitalDistributionMain, 'file')==0
    mkdir(DirectoryIntitalDistributionMain); end
if exist(DirectoryIntitalDistribution, 'file')==0
    mkdir(DirectoryIntitalDistribution); end

if MicrohabitatType==1
    DirectoryMicrohabitat=strcat(DirectoryMicrohabitatMain,'\DynamicForests\',FolderMicrohabitat,...
        '\Microhabitat_',FolderMicrohabitat,'_Rep',num2str(Replicate));
elseif MicrohabitatType==2
    DirectoryMicrohabitat=strcat(DirectoryMicrohabitatMain,'\StaticForests\',FolderMicrohabitat,...
        '\Microhabitat_',FolderMicrohabitat,'_Rep',num2str(Replicate));
elseif MicrohabitatType==3
    DirectoryMicrohabitat=strcat(DirectoryMicrohabitatMain,'\UniformForests\',FolderMicrohabitat,...
        '\Microhabitat_',FolderMicrohabitat,'_Rep',num2str(Replicate));
end

DirectorySpeciesPools=strcat(DirectorySpeciesPoolsMain,'\',FolderSpeciesPools);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Load parameters saved along with the microhabitat and species pool files

%Load plot dimensions if an artifical theoretical forest is used
load(strcat(DirectoryMicrohabitat,'\dimPlot.mat'))

%Calculate individuals per species if normalization per hectare (ScalingPerHa=1) is chosen
if ScalingPerHa==1
    IndividualsPerSpecies=IndividualsPerSpecies*((dimPlot(1)*dimPlot(2))/10000);
end

%Get number of species from species pool file
Input_file=strcat(DirectorySpeciesPools,'\SpeciesPool',num2str(numSpeciesPools(1)),'.csv');
SpeciesPool=dlmread(Input_file,'\t');
NumberSpecies=size(SpeciesPool,1);

%Get number of total individuals for eacg replicate
TotalIndividuals=NumberSpecies*IndividualsPerSpecies;
NumberMaturesPerSpecies=round(IndividualsPerSpecies*(PercentageMaturePerSpecies/100));

%Load initial microhabitat matrix
FileInitalMatrix=strcat(DirectoryMicrohabitat,'\MicrohabitatMatrix',num2str(TimeStep),'.mat');
load(FileInitalMatrix)

%Set real light values (in microhabitat, relative light values are saved)
Microhabitat(:,:,:,3)=Microhabitat(:,:,:,3).*Imax;

[Test, ColumnHeaders]=xlsread(strcat(DirectorySpeciesPools,'\ColumnHeaders.xls'));
ColumnHeaders={ColumnHeaders{:},'X','Y','Z','Mass','Status','IndividualID','SurfaceAreaOccupied','Age'};
xlswrite(strcat(DirectoryIntitalDistribution,'\ColumnHeaders.xls'),ColumnHeaders)

%Get numbers of columns used in this script
ColMinLight=find(ismember(ColumnHeaders,'MinLight'));
ColMaxLight=find(ismember(ColumnHeaders,'MaxLight'));

%Save directory of the species pool in seperate file
dlmwrite(strcat(DirectoryIntitalDistribution,'\DirectorySpeciesPool.txt'),DirectorySpeciesPools,''); %Save associated directory of species pools

%Save directory of the microhabitat matrix in seperate file
dlmwrite(strcat(DirectoryEpiphyteModel,'\DirectoryMicrohabitatMatrix.txt'),DirectoryMicrohabitat,''); %Save associated directory of species pools

%Growth function. Used here to approximate the age of the individuals
MassFunctionOfAge=@(MaxMass,K,Age) (MaxMass*(1-exp(-K*(Age))));
AgeFunctionOfMass=@(MaxMass,Mass,K) (-log(1-(Mass/MaxMass))/K);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Main loop for Single Species Model
if SingleSpeciesModel==1
    
    AvailableSurfaceAreaForSpecies=zeros(dimPlot(1),dimPlot(2),dimPlot(3),NumberSpecies,'single');  
    PotentialVoxelsForSpecies=zeros(dimPlot(1)*dimPlot(2)*dimPlot(3),3,NumberSpecies,'int8');
    PotentialVoxelsForIndividual=zeros(dimPlot(1)*dimPlot(2)*dimPlot(3),3,'int8');

    for numPool=numSpeciesPools(1):numSpeciesPools(2)

        fprintf('Number species pool: %d \n', numPool);
        %Laden des species pools
        Input_file=strcat(DirectorySpeciesPools,'\','SpeciesPool',num2str(numPool),'.csv');
        SpeciesPool=dlmread(Input_file,'\t');

        SizeSpeciesPool=size(SpeciesPool);
        

        for numReplicates=1:replicatePerSpeciesPool
            
            fprintf('Number replicate: %d \n', numReplicates);
            %Initialize epiphyte matrix
            IntitalEpiphyteMatrix=zeros(TotalIndividuals,(size(SpeciesPool,2)+7));
            
            %Initialize Available surfavce area
            AvailableSurfaceArea=Microhabitat(:,:,:,1); %Matrix to trace the still available surface area per voxel

            %Fill InitialEpiphyteMatrix with species trait informations and the
            %initial size of each individual
            for numSpecies=1:NumberSpecies
                for numIndividual=1:IndividualsPerSpecies

                    %Copy trait data from SpeciesPool to InitalEpiphyteMatrix
                    IntitalEpiphyteMatrix(((numSpecies-1)*IndividualsPerSpecies)+numIndividual,1:SizeSpeciesPool(2))=SpeciesPool(numSpecies,:);

                    %Get size of individual
                    if numIndividual<=NumberMaturesPerSpecies
                        SizeOfIndividual=random('unif',SpeciesPool(numSpecies,3),SpeciesPool(numSpecies,2));%Size of mature individuals
                    else
                        SizeOfIndividual=random('unif',0,SpeciesPool(numSpecies,3));%Size of juvenile individuals
                    end

                    %Store initial size of individual
                    IntitalEpiphyteMatrix(((numSpecies-1)*IndividualsPerSpecies)+numIndividual,SizeSpeciesPool(2)+4)=...
                        SizeOfIndividual;
                    
                    %store initial age of individual (age when it would have grown under optimal conditions)
                    IntitalEpiphyteMatrix(((numSpecies-1)*IndividualsPerSpecies)+numIndividual,SizeSpeciesPool(2)+8)=...
                        round(AgeFunctionOfMass(SpeciesPool(numSpecies,2),SizeOfIndividual,SpeciesPool(numSpecies,4)));

                end
            end
            
            %Store individual ID for each individual
            IntitalEpiphyteMatrix(1:TotalIndividuals,(SizeSpeciesPool(2)+6))=1:TotalIndividuals;

            %Calculate the surface area needed to support an individual of
            %this size =SurfaceAreaNeededInVoxel
            IntitalEpiphyteMatrix(:,(SizeSpeciesPool(2)+7))=(IntitalEpiphyteMatrix(:,(SizeSpeciesPool(2)+4)).^(2/3))/SurfaceBiomassScaling;


            %loop over all species
            for numSpecies=1:NumberSpecies
            
                %Get subset of indiduals for each species
                IntitalEpiphyteMatrixSub=IntitalEpiphyteMatrix(IntitalEpiphyteMatrix(:,1)==numSpecies,:);
                
                %loop through all individuals, beginning with the largest (competition)
                %sort IntitalEpiphyteMatrixSub by size
                IntitalEpiphyteMatrixSub=sortrows(IntitalEpiphyteMatrixSub,-(SizeSpeciesPool(2)+4));
                NumNoSurface=0;

                %Calculate potential voxels for for each species which fullfil
                %their niche requirments (to save time they are precomputed here)
                [x,y,z]=ind2sub(size(Microhabitat),find(Microhabitat(:,:,:,1)>0 ...
                 &  Microhabitat(:,:,:,3)>=SpeciesPool(numSpecies,ColMinLight)  & Microhabitat(:,:,:,3)<=SpeciesPool(numSpecies,ColMaxLight)));

                for i=1:length(IntitalEpiphyteMatrixSub)

                    randNumbers=randperm(length(x));

                    for PotVoxels=1:length(x)

                        if(AvailableSurfaceArea(x(randNumbers(PotVoxels)),y(randNumbers(PotVoxels)),z(randNumbers(PotVoxels)))...
                                            >IntitalEpiphyteMatrixSub(i,(SizeSpeciesPool(2)+7)))

                            IntitalEpiphyteMatrixSub(i,SizeSpeciesPool(2)+1)=x(randNumbers(PotVoxels));
                            IntitalEpiphyteMatrixSub(i,SizeSpeciesPool(2)+2)=y(randNumbers(PotVoxels));
                            IntitalEpiphyteMatrixSub(i,SizeSpeciesPool(2)+3)=z(randNumbers(PotVoxels));

                            %Set status of individual: status=1 => alive
                            IntitalEpiphyteMatrixSub(i,SizeSpeciesPool(2)+5)=1;

                            AvailableSurfaceArea(x(randNumbers(PotVoxels)),y(randNumbers(PotVoxels)),z(randNumbers(PotVoxels)))=...
                                AvailableSurfaceArea(x(randNumbers(PotVoxels)),y(randNumbers(PotVoxels)),z(randNumbers(PotVoxels)))-...
                                (IntitalEpiphyteMatrixSub(i,(SizeSpeciesPool(2)+7)));         

                            break;
                        end
                    end
                end

                %Set the coordinates to 1 for all individuals that did not find any
                %suitable habitat
                IntitalEpiphyteMatrixSub(IntitalEpiphyteMatrixSub(:,SizeSpeciesPool(2)+1)==0,SizeSpeciesPool(2)+5)=2;
                IntitalEpiphyteMatrixSub(IntitalEpiphyteMatrixSub(:,SizeSpeciesPool(2)+1)==0,SizeSpeciesPool(2)+1)=1;
                IntitalEpiphyteMatrixSub(IntitalEpiphyteMatrixSub(:,SizeSpeciesPool(2)+2)==0,SizeSpeciesPool(2)+2)=1;
                IntitalEpiphyteMatrixSub(IntitalEpiphyteMatrixSub(:,SizeSpeciesPool(2)+3)==0,SizeSpeciesPool(2)+3)=1;
                
                IntitalEpiphyteMatrix(((numSpecies-1)*IndividualsPerSpecies)+1:(numSpecies*IndividualsPerSpecies),:)=IntitalEpiphyteMatrixSub;
                
            end

            %Save Inital Epiphyte Matrix
            SaveFile=strcat(DirectoryIntitalDistribution,'\ID_SpeciesP_',num2str(numPool),'_Rep_',num2str(numReplicates),'.csv');
            mkdir(DirectoryIntitalDistribution)
            dlmwrite(SaveFile,IntitalEpiphyteMatrix,'\t')

        end

    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Main loop for community model
if SingleSpeciesModel==0
    for numPool=numSpeciesPools(1):numSpeciesPools(2)

        fprintf('Number species pool: %d \n', numPool);

        %Laden des species pools
        Input_file=strcat(DirectorySpeciesPools,'\','SpeciesPool',num2str(numPool),'.csv');
        SpeciesPool=dlmread(Input_file,'\t');
        SizeSpeciesPool=size(SpeciesPool);

        for numReplicates=1:replicatePerSpeciesPool

            fprintf('Number replicate: %d \n', numReplicates);
            %Initialize epiphyte matrix
            IntitalEpiphyteMatrix=zeros(TotalIndividuals,(size(SpeciesPool,2)+7));
            
            %Initialize Available surfavce area
            AvailableSurfaceArea=Microhabitat(:,:,:,1); %Matrix to trace the still available surface area per voxel

            %Fill InitialEpiphyteMatrix with species trait informations and the
            %initial size of each individual
            for numSpecies=1:NumberSpecies
                for numIndividual=1:IndividualsPerSpecies

                    %Copy trait data from SpeciesPool to InitalEpiphyteMatrix
                    IntitalEpiphyteMatrix(((numSpecies-1)*IndividualsPerSpecies)+numIndividual,1:SizeSpeciesPool(2))=SpeciesPool(numSpecies,:);

                    %Get size of individual
                    if numIndividual<=NumberMaturesPerSpecies
                        SizeOfIndividual=random('unif',SpeciesPool(numSpecies,3),SpeciesPool(numSpecies,2));%Size of mature individuals
                    else
                        SizeOfIndividual=random('unif',0,SpeciesPool(numSpecies,3));%Size of juvenile individuals
                    end

                    %Store initial size of individual
                    IntitalEpiphyteMatrix(((numSpecies-1)*IndividualsPerSpecies)+numIndividual,SizeSpeciesPool(2)+4)=...
                        SizeOfIndividual;
                    
                    %store initial age of individual (age when it would have grown under optimal conditions)
                    IntitalEpiphyteMatrix(((numSpecies-1)*IndividualsPerSpecies)+numIndividual,SizeSpeciesPool(2)+8)=...
                        round(AgeFunctionOfMass(SpeciesPool(numSpecies,2),SizeOfIndividual,SpeciesPool(numSpecies,4)));

                end
            end
            
            %Store individual ID for each individual
            IntitalEpiphyteMatrix(1:TotalIndividuals,(SizeSpeciesPool(2)+6))=1:TotalIndividuals;

            %Calculate the surface area needed to support an individual of
            %this size =SurfaceAreaNeededInVoxel
            IntitalEpiphyteMatrix(:,(SizeSpeciesPool(2)+7))=(IntitalEpiphyteMatrix(:,(SizeSpeciesPool(2)+4)).^(2/3))/SurfaceBiomassScaling;

            %loop randomly through all individuals and select suitable
            %voxel for each. The suitable voxel can either be the voxel
            %with the highest available surface area (MethodVoxel=1), or a random voxel
            %(MethodVoxel=0)
            MethodVoxel=0;
            RandNumInd=randsample(TotalIndividuals,TotalIndividuals);
             
            for i=1:TotalIndividuals
               
                disp(['Individual number ', num2str(i)]);
                NumIndRand=RandNumInd(i);
                
                %Find all suitable voxels for this individual
                MinLightInd=IntitalEpiphyteMatrix(NumIndRand,ColMinLight);
                MaxLightInd=IntitalEpiphyteMatrix(NumIndRand,ColMaxLight);
                AreaNeededInd=IntitalEpiphyteMatrix(NumIndRand,(SizeSpeciesPool(2)+7));
             
                %1. Get the postions of all voxels fullfilling the
                %requirements of the individual (light+area)
                SuitableVoxels=find(AvailableSurfaceArea(:,:,:)>AreaNeededInd ...
                    &  Microhabitat(:,:,:,3)>=MinLightInd  & Microhabitat(:,:,:,3)<=MaxLightInd);  

                 %Choose one of the suitable voxels based on the specified
                %Method (if suitable voxels are available)
                if size(SuitableVoxels,1)>0
                    if MethodVoxel==1
                        MaxVal=find(AvailableSurfaceArea(SuitableVoxels)==max(AvailableSurfaceArea(SuitableVoxels)));
                        [x,y,z]=ind2sub(size(AvailableSurfaceArea),SuitableVoxels(MaxVal(1)));
                    elseif MethodVoxel==0
                        RandVal=randi(size(SuitableVoxels,1));
                        [x,y,z]=ind2sub(size(AvailableSurfaceArea),SuitableVoxels(RandVal));
                    end

                    %Update available Surface Area
                    AvailableSurfaceArea(x,y,z)=AvailableSurfaceArea(x,y,z)-AreaNeededInd;

                    %Update Inital Epiphyte Matrix
                    IntitalEpiphyteMatrix(NumIndRand,SizeSpeciesPool(2)+1)=x;
                    IntitalEpiphyteMatrix(NumIndRand,SizeSpeciesPool(2)+2)=y;
                    IntitalEpiphyteMatrix(NumIndRand,SizeSpeciesPool(2)+3)=z;

                    %Set status of individual: status=1 => alive
                    IntitalEpiphyteMatrix(NumIndRand,SizeSpeciesPool(2)+5)=1;
                
                else
                    
                    %Set status of individual: status=2 =>> dead
                    IntitalEpiphyteMatrix(NumIndRand,SizeSpeciesPool(2)+5)=2;
                    
                    %Set coordinates to 1 (might cause problems in later model if not)
                    IntitalEpiphyteMatrix(NumIndRand,SizeSpeciesPool(2)+1)=1;
                    IntitalEpiphyteMatrix(NumIndRand,SizeSpeciesPool(2)+2)=1;
                    IntitalEpiphyteMatrix(NumIndRand,SizeSpeciesPool(2)+3)=1;
                end
            end
            
            %Summary
            IndividualsWithoutVoxels=size(find(IntitalEpiphyteMatrix(:,SizeSpeciesPool(2)+5)==2),1)
            
            %Save Inital Epiphyte Matrix
            SaveFile=strcat(DirectoryIntitalDistribution,'\ID_SpeciesP_',num2str(numPool),'_Rep_',num2str(numReplicates),'.csv');
            dlmwrite(SaveFile,IntitalEpiphyteMatrix,'\t')
        end
    end
end
