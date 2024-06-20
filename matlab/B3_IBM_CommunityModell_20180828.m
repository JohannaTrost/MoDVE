%Epiphte IBM - Model
%This model simulates the development of the entire epiphyte community
clear all
close all
clc

%turn warning messages off
warning('off','all')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Parameters that need to be specified/checked before running this script

%Folder of epiphyte models (these models are simulated in this order)
%The names of the models in the Folder "EpiphyteModels" are needed here
%FolderEpiphyteModels={'Dynamic_ForestModel_EM_20150914_2'};
FolderEpiphyteModels={'ForestModel_LoggedG35_50x50_Rep4'};
MicrohabitatType=1; %Define which type of forest the microhabitat belongs to. 1: dynamic forest, 2: static forest, 3: uniform forest
SingleSpeciesModel=0; %1: Single species model, 0: Community model

%Choose initial distributions (have to be located in 'FolderEpiphyteModel\IniDist\')
FolderInitialDistributions={'SP_Random_IA_2_IR_60_TimeS_200_Pop1010-1020'};

DirectoryModelMain='C:\Gunnar\EpiphyteModels\EpiphyteModels';

%Model parameters
timeSteps=600; %Model for timeSteps beginning at the time step given by the initial distribution

%Density of individuals per ha at which to stop the simulationof the community and
%move to the next replicate (to prevent exploding communities)
StopCriterionHa=3000000; %Individuals per ha

%Choose species pools to use and number of replicates per species pool
numSpeciesPools=[1,10]; %Start and end number of  species pools (if the species pools do not exist, they are automatically skipped)replicatePerSpeciesPool=5; %Number of replicates per species pool  (if the replicates do not exist, they are automatically skipped)
replicatePerSpeciesPool=1; %Number of replicates per species pool  (if the replicates do not exist, they are automatically skipped)

SurfaceBiomassScaling=10000*0.01; %cm^2 per m^2
Imax=900; %maximum light above canopy

%Competition Methods; defines which individuals are removed in voxels which
%are entirely filled. 1:size (small individuals are outcompetet by larger ones); 2:random competition
CompetitionMethod=1;

%Mortality method (complete random or scaling with mass according to metabolic theory);
MortalityMethod=1; %0: random mortality; 1: scaling with mass to the exponent -1/4
MortRateRandom=0.1;
MortRateMass=0.1;
MortRateMassScaling=-0.25; %widely used scaling fator
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for MicrohabitatNumber=1:length(FolderEpiphyteModels);

    for InititalDistNumber=1:length(FolderInitialDistributions);

        %Epiphyte model
        FolderEpiphyteModel=char(FolderEpiphyteModels{MicrohabitatNumber});

        %Initial distribution
        FolderInitialDistribution=char(FolderInitialDistributions{InititalDistNumber});

        %Get initial time step
        FolderInitialDistributionsTemp=char(FolderInitialDistributions(InititalDistNumber));
        PosUnderscores=regexp(FolderInitialDistributionsTemp, '_');
        InitialTimeStepTemp=FolderInitialDistributionsTemp((PosUnderscores(7)+1):(PosUnderscores(8)-1));
        InitialTimeStep=str2double(InitialTimeStepTemp);

        if MicrohabitatType==1
           FolderEpiphyteModelMain='DynamicForests';
        elseif MicrohabitatType==2
           FolderEpiphyteModelMain='StaticForests';
        elseif MicrohabitatType==3
           FolderEpiphyteModelMain='UniformForests';
        end

        if SingleSpeciesModel==1
            FolderModelType='SingleSpeciesModels';
        elseif SingleSpeciesModel==0
            FolderModelType='CommunityModels';
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Get main model folder and folder where initial distribution is saved
        DirectoryEpiphyteModel=strcat(DirectoryModelMain,'\',FolderModelType,'\',FolderEpiphyteModelMain,'\',FolderEpiphyteModel);
        DirectoryIntitalDistribution=strcat(DirectoryEpiphyteModel,'\IniDist\',FolderInitialDistribution);

        %Get folder where microhabitat and species pool is saved
        fid=fopen(strcat(DirectoryEpiphyteModel,'\DirectoryMicrohabitatMatrix.txt')); %The directory of the micorhabitat matrix is stored here
        MicMat = textscan(fid, '%s','Delimiter','\t');
        DirectoryMicrohabitat=char(MicMat{1});
        fclose(fid);

        fid=fopen(strcat(DirectoryIntitalDistribution,'\DirectorySpeciesPool.txt')); %The directory of the species pool is stored here
        SpePool = textscan(fid, '%s','Delimiter','\t');
        DirectorySpeciesPools=char(SpePool{1});
        fclose(fid);

        %Create folder to save the model results
        DirectoryModelResults=strcat(DirectoryEpiphyteModel,'\ModelResults\',FolderInitialDistribution);
        mkdir(DirectoryModelResults)

        %Load plot dimensions
        load(strcat(DirectoryMicrohabitat,'\dimPlot.mat'))

        %Set StopCriterion for this simulation
        StopCriterion=StopCriterionHa*(dimPlot(1)*dimPlot(2)*0.0001);

        %Load TraitRanges (ranges used to create the species pool)
        FileTraitRanges=strcat(DirectorySpeciesPools,'\TraitRanges.csv');
        TraitRanges=dlmread(FileTraitRanges, '\t');

        %Get information if random, sequential or neutral  species pool was created (this affects the way how to
        %get the information from the trait ranges (see below)
        load(strcat(DirectorySpeciesPools,'\SpeciesPoolType.mat'));

        if  SpeciesPoolType==0 %random species pool
            SlopeRecruitment=TraitRanges(1,1);
            InterceptRecruitment=TraitRanges(2,1);
        elseif SpeciesPoolType==1 %sequential species pool
            SlopeRecruitment=TraitRanges(1,1);
            InterceptRecruitment=TraitRanges(2,1);
        elseif SpeciesPoolType==2 %neutral species pool
            SlopeRecruitment=TraitRanges(1,1);
            InterceptRecruitment=TraitRanges(2,1);
        end

          %Load column header and assign them accordingly
        [Test, ColumnHeaders]=xlsread(strcat(DirectoryIntitalDistribution,'\ColumnHeaders.xls'));
        [Test, ColumnHeadersSpeciesPool]=xlsread(strcat(DirectorySpeciesPools,'\ColumnHeaders.xls'));
        NumColSpeciesPool=length(ColumnHeadersSpeciesPool);

        ColSpeciesID=find(ismember(ColumnHeaders,'SpeciesID'));
        ColMaximumMass=find(ismember(ColumnHeaders,'MaximumMass'));
        ColMassAtMaturity=find(ismember(ColumnHeaders,'MassAtMaturity'));
        ColGrowthRate=find(ismember(ColumnHeaders,'GrowthRate'));
        ColDispersalKernel=find(ismember(ColumnHeaders,'DispersalKernel'));
        ColDispersalKernelAsymmetry=find(ismember(ColumnHeaders,'DispersalKernelAsymmetry'));
        ColRecruitmentInvestmentRel=find(ismember(ColumnHeaders,'RecruitmentInvestmentRel'));
        ColRecruitmentInc=find(ismember(ColumnHeaders,'RecruitmentInc'));
        ColMinLight=find(ismember(ColumnHeaders,'MinLight'));
        ColMaxLight=find(ismember(ColumnHeaders,'MaxLight'));
        ColOptimumLight=find(ismember(ColumnHeaders,'OptimumLight'));
        ColLightBreadth=find(ismember(ColumnHeaders,'LightBreadth'));
        ColLightResponseA=find(ismember(ColumnHeaders,'LightResponseA'));
        ColLightResponseB=find(ismember(ColumnHeaders,'LightResponseB'));
        ColLightResponseC=find(ismember(ColumnHeaders,'LightResponseC'));
        ColX=find(ismember(ColumnHeaders,'X'));
        ColY=find(ismember(ColumnHeaders,'Y'));
        ColZ=find(ismember(ColumnHeaders,'Z'));
        ColMass=find(ismember(ColumnHeaders,'Mass'));
        ColStatus=find(ismember(ColumnHeaders,'Status'));
        ColIndividualID=find(ismember(ColumnHeaders,'IndividualID'));
        ColSurfaceAreaOccupied=find(ismember(ColumnHeaders,'SurfaceAreaOccupied'));
        ColAge=find(ismember(ColumnHeaders,'Age'));
        ColTotalSurfaceInVoxel=length(ColumnHeaders)+1;
        ColLightInVoxel=length(ColumnHeaders)+2;
        ColSurfaceLossInVoxel=length(ColumnHeaders)+3;
        TotalColsEpiphyteMatrix=ColSurfaceLossInVoxel;

        %The following information are saved for each individual (the traits can be accessed via the SpeciesID)
        ColumsToSave=[ColSpeciesID,ColIndividualID,ColStatus,ColMass,ColAge,ColX,ColY,ColZ,...
            ColTotalSurfaceInVoxel,ColSurfaceLossInVoxel,ColLightInVoxel];
        %Headers of this matrix
        SummaryMatrixIndividualsHeaders={'SpeciesID','IndividualID','Status','Mass','Age',...
            'X','Y','Z','TotalSurfaceInVoxel','SurfaceLossInVoxel','LightInVoxel'};

        if length(SummaryMatrixIndividualsHeaders)~=length(ColumsToSave)
            disp('Headers of individual matrix do not match with number of colums => rework!!!!')
            return
        end

        %Following the columns in the species matrix refering to a trait or
        %variable. This is handy if the epiphyte matrix changes
        ColSSpeciesID=1;
        ColSNumberIndividualsBeginning=2;
        ColSNumberIndividualsEnd=3;
        ColSNumberMatureIndividuals=4;
        ColSNumberRecruits=5;
        ColSNumberRecruitsPotential=6;
        ColSNumberMortalityBranchFall=7;
        ColSNumberMortalityLight=8;
        ColSNumberMortalityCompetition=9;
        ColSNumberMortalityNatural=10;
        ColSNumberPopulationGrowthRate=11;
        ColSNumberPopulationGrowthRateLog=12;
        ColSNumberBirthRate=13;
        ColSNumberDeathRate=14;
        ColSAverageSize=15;
        ColSAverageAge=16;
        ColSMinLight=17;
        ColSMaxLight=18;
        ColSMeanLight=19;
        ColSMinHeight=20;
        ColSMaxHeight=21;
        ColSMeanHeight=22;
        TotalColsSpeciesMatrix=ColSMeanHeight;

        %Headers of matrix
        SummaryMatrixSpeciesHeaders={'TimeStep','SpeciesID','NumberIndividualsBeginning','NumberIndividualsEnd','NumberMatureIndividuals',...
            'NumberRecruits','NumberRecruitsPotential','NumberMortalityBranchFall','NumberMortalityLight',...
            'NumberMortalityCompetition','NumberMortalityNatural',...
            'PopulationGrowthRate','PopulationGrowthRateLog','BirthRate','DeathRate',...
            'AverageMass','AverageAge','MinLight','MaxLight','MeanLight',...
            'MinHeight','MaxHeight','MeanHeight'};

        if length(SummaryMatrixSpeciesHeaders)~=(TotalColsSpeciesMatrix+1)
            disp('Headerd of species matrix do not match with number of colums => rework!!!!')
            return
        end

        %Headers of matrix
        SummaryMatrixCommunityHeaders={'timeStep','NumberSpeciesBeginning','NumberSpeciesEnd','NumberIndividualsBeginning',...
            'NumberIndividualsEnd','Recruits','MortalityBranchFall','MortalityLight',...
            'MortalityCompetition','MortalityNatural','BranchSurfaceIndex','EpiphyteFilling'};

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Functions used in the model

        %Bertalanffy Growth
        GrowthRate=@(MaxMass,Mass,K) (K*(MaxMass-Mass));

        %Parabolic Optimum function
        Parabol=@(a,b,c,x) a*x^2+b*x+c;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Main loop for the community model

        SummaryMatrixIndividualsHeaders={'TimeStep' SummaryMatrixIndividualsHeaders{:}};

        for numPool=numSpeciesPools(1):numSpeciesPools(2)

            %Check if a initial distribution for the species pool exists. If not, move on to the next species pool
            FileNameInitalDistributionPool=strcat(DirectoryIntitalDistribution,'\ID_SpeciesP_',num2str(numPool),'_Rep_1.csv');
            if exist(FileNameInitalDistributionPool, 'file')==0
                break;
            end

            %First step: create probability matrices for each species
            %Load species pool
            SpeciesPool=dlmread(strcat(DirectorySpeciesPools,'\SpeciesPool',num2str(numPool),'.csv'),'\t');
            NumberOfSpecies=size(SpeciesPool,1); %number of species per 25X25m plot

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Erzeugen der Distanzmatrix und der Wahrscheinlichkeitsmatrix f�r jede Art
            %Dimensionen der Dispersal matrix
            dimX=dimPlot(1)*2+1;
            dimY=dimPlot(2)*2+1;
            dimZ=dimPlot(3)*2+1;

            %Erzeugen der Distanzmatrix mit allen Distanzen zum
            centralPoint=[floor(dimX/2)+1, floor(dimY/2)+1, floor(dimZ/2)+1 ];
            DistanceMatrix=zeros(dimX,dimY,dimZ);

            for i=1:dimX
                for j=1:dimY
                    for k=1:dimZ
                        DistanceMatrix(i,j,k)=pdist([i,j,k; centralPoint(1),centralPoint(2),centralPoint(3)]);
                    end
                end
            end

            %Erzeugen der Wahrschienlichkeitsmatrix anhand der Distanzmatrix und dem
            %artspezischien Wert bb aus der Epiphytenmatrix
            negExp = @(distance,bb) exp(-distance.*bb); %Negative Exponential function
            ProbabilityMatrix=zeros(dimX,dimY,dimZ,NumberOfSpecies);
            ProbabilityMatrixNormalized=zeros(dimX,dimY,dimZ,NumberOfSpecies,'double');

            for i=1:NumberOfSpecies
                exponentE=SpeciesPool(i,ColDispersalKernel);
                dispersalAsymmetry=SpeciesPool(i,ColDispersalKernelAsymmetry);

                ProbabilityMatrix(:,:,:,i)=negExp(DistanceMatrix(:,:,:),exponentE);

                %Apply dispersal asymmetry (probability to disperse downwards higher than upwards dispersal)
                ProbabilityMatrix(:,:,centralPoint(3):dimZ)=ProbabilityMatrix(:,:,centralPoint(3):dimZ).*((1-dispersalAsymmetry)/0.5);
                ProbabilityMatrix(:,:,1:(centralPoint(3)-1))=ProbabilityMatrix(:,:,1:(centralPoint(3)-1)).*(dispersalAsymmetry/0.5);

                ProbabilityMatrixNormalized(:,:,:,i)=ProbabilityMatrix(:,:,:,i)./sum(sum(sum(ProbabilityMatrix(:,:,:,i))));
            end
            clear ProbabilityMatrix %to save memory

            %Main model loop for each replicate
            for r=1:replicatePerSpeciesPool

                %Check if a initial distribution for the species pool exists. If not, move on to the next species pool
                FileNameInitalDistribution=strcat(DirectoryIntitalDistribution,'\ID_SpeciesP_',num2str(numPool),'_Rep_',num2str(r),'.csv');
                if exist(FileNameInitalDistribution, 'file')==0
                    break;
                end

                %Create Save-Directory for each each replicate/initialDistribution
                DirectoryModelResultsRun=strcat(DirectoryModelResults,'\ID_SpeciesP_',num2str(numPool),'_Rep_',num2str(r));
                mkdir(DirectoryModelResultsRun)

                %Load initial epiphyte distribution
                E=dlmread(strcat(DirectoryIntitalDistribution,'\ID_SpeciesP_',num2str(numPool),'_Rep_',num2str(r),'.csv'),'\t'); %E for epiphytes

                %Add column to E for additional information
                E=[E zeros(size(E,1),4)];

                MaxIndividualID=size(E,1); %to trace individual IDs

                %Initialize Matrix where community parameters are save
                SummaryMatrixCommunity=zeros(timeSteps,length(SummaryMatrixCommunityHeaders));

                %Load microhabitat matrix if a uniform or static forest is simulated (only needs to be loaded once an not envery timestep)
                if MicrohabitatType==2 || MicrohabitatType==3
                    load(strcat(DirectoryMicrohabitat,'\MicrohabitatMatrix1.mat'))
                    Microhabitat(:,:,:,3)=Microhabitat(:,:,:,3)*Imax; %In the microhabitat matrix, the realtive light extinction is stored: convert to light values in ?mol*m-2*s-1
                    pot_habitat=zeros(size(Microhabitat,1),size(Microhabitat,2),size(Microhabitat,3));
                end

                %Initialize matrices where the aggregated information on species level are saved
                SummaryMatrixSpeciesSave=zeros(timeSteps*NumberOfSpecies,TotalColsSpeciesMatrix);

                %Save directory where species pool is stored
                dlmwrite(strcat(DirectoryModelResultsRun,'\DirectorySpeciesPool.txt'),DirectorySpeciesPools,''); %Save associated directory of species pools


                %Save Column header for further use
                fid = fopen(strcat(DirectoryModelResultsRun,'\HeadersCommunitySummary.csv'),'w');
                fprintf(fid,[ repmat('%s,',1,size(SummaryMatrixCommunityHeaders,2)-1),'%s\n'],SummaryMatrixCommunityHeaders{1,:});
                fclose(fid);

                fid = fopen(strcat(DirectoryModelResultsRun,'\HeadersSpeciesSummary.csv'),'w');
                fprintf(fid,[ repmat('%s,',1,size(SummaryMatrixSpeciesHeaders,2)-1),'%s\n'],SummaryMatrixSpeciesHeaders{1,:});
                fclose(fid);

                fid = fopen(strcat(DirectoryModelResultsRun,'\HeadersIndividualMatrix.csv'),'w');
                fprintf(fid,[ repmat('%s,',1,size(SummaryMatrixIndividualsHeaders,2)-1),'%s\n'],SummaryMatrixIndividualsHeaders{1,:});
                fclose(fid);

                %Initialize Matrix where speceies parameters are save
                SummaryMatrixSpecies=zeros(timeSteps*NumberOfSpecies,TotalColsSpeciesMatrix);

                for t=InitialTimeStep:(InitialTimeStep+timeSteps)

                   %Check if the stop criterion is met
                   if size(E(E(:,ColStatus)==1),1)>StopCriterion
                        break;
                   end

                   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %1. Dispersal

                    %Load microhabitat matrix for specific timeStep if dynamic forest is simulated
                    if MicrohabitatType==1
                        load(strcat(DirectoryMicrohabitat,'\MicrohabitatMatrix',num2str(t),'.mat'))
                        Microhabitat(:,:,:,3)=Microhabitat(:,:,:,3)*Imax; %In the microhabitat matrix, the realtive light extinction is stored: convert to light values in ?mol*m-2*s-1
                        pot_habitat=zeros(size(Microhabitat,1),size(Microhabitat,2),size(Microhabitat,3));
                    end

                    %Store number of individuals at beginning of time step
                    for g=1:NumberOfSpecies
                        IntialNumberIndividuals(g)=size(E(E(:,ColSpeciesID)==g & E(:,ColStatus)==1 ,:),1);
                    end
                    IntialNumberIndividualsTotal=size(E(E(:,ColStatus)==1 ,:),1);
                    InitialNumberSpecies=length(unique(E(E(:,ColStatus)==1,ColSpeciesID)));
                    NumberRecruitsPerSpecies=zeros(1,NumberOfSpecies);

                    %Calculate free surface area per voxel
                    AvailableSurfaceArea=Microhabitat(:,:,:,1);
                    for i=1:size(E,1)
                        SurfaceAreaNeededInVoxel=(E(i,ColMass)^(2/3))/SurfaceBiomassScaling;
                        AvailableSurfaceArea(E(i,ColX),E(i,ColY),E(i,ColZ))=max(0,(AvailableSurfaceArea(E(i,ColX),E(i,ColY),E(i,ColZ))...
                            -SurfaceAreaNeededInVoxel));
                    end

                    if ~isempty(E) %Check if there are species left

                        unique_species=unique(E(:,1)); %list with species IDs of all present species

                        %loop over all species
                        for i=1:length(unique_species) %Number of species

                            %Generate initially empty matrix to store the probabilities for recruitment
                            ProbabilityMatrixPerSpecies=zeros(dimPlot(1),dimPlot(2),dimPlot(3));

                            %Matrix containing all mature individuals of one species
                            MatureIndividulsPerSpecies=E(E(:,ColSpeciesID)==unique_species(i) & E(:,ColMass)>=E(:,ColMassAtMaturity) ,1:ColAge);

                            if ~isempty(MatureIndividulsPerSpecies)

                                %Probability matrix for each species: Depending on the position of each mature individual, the total
                                %probability for the species is calculated.
                                %The second part of the equation accounts for the actual size of the individual
                                %in relation to the maximum size for which the recruitment per individual is defined
                                for j=1:size(MatureIndividulsPerSpecies,1) %Number of mature individuals per species
                                    ProbabilityMatrixPerSpecies=ProbabilityMatrixPerSpecies+...
                                    ProbabilityMatrixNormalized(centralPoint(1)-MatureIndividulsPerSpecies(j,ColX)+...
                                    1:centralPoint(1)-MatureIndividulsPerSpecies(j,ColX)+dimPlot(1),...
                                    centralPoint(2)-MatureIndividulsPerSpecies(j,ColY)+...
                                    1:centralPoint(2)-MatureIndividulsPerSpecies(j,ColY)+dimPlot(2),...
                                    centralPoint(3)-MatureIndividulsPerSpecies(j,ColZ)+...
                                    1:centralPoint(3)-MatureIndividulsPerSpecies(j,ColZ)+dimPlot(3),MatureIndividulsPerSpecies(j,ColSpeciesID)).*...
                                    ((InterceptRecruitment)*...
                                    MatureIndividulsPerSpecies(j,ColRecruitmentInvestmentRel))
                                    ;
                                %Here I change this part to the simple
                                %version :needs to be tested!!!

                                end

                                %Store potential normalized number of recruits in SummaryMatrixSpecies
                                SummaryMatrixSpecies(((i-1)*timeSteps)+t,ColSNumberRecruitsPotential)=...
                                    sum(sum(sum(ProbabilityMatrixPerSpecies))); %potential recruitment

                                %Matix containing all voxel for which the light requirements are fulfilled
                                pot_habitat(:,:,:)=Microhabitat(:,:,:,3)>=...
                                   MatureIndividulsPerSpecies(1,ColMinLight) &...
                                   Microhabitat(:,:,:,3)<=MatureIndividulsPerSpecies(1,ColMaxLight);

                                %Final probabiliy matrix for new recruits
                                probability_recruits=ProbabilityMatrixPerSpecies...
                                    .*pot_habitat(:,:,:)...
                                    .*AvailableSurfaceArea(:,:,:);

                                %Calculate number of recuits based on final probability matrix
                                Recruits=poissrnd(probability_recruits);

                                %Add new recruits to epiphyte matrix
                                num_recruits=sum(sum(sum(Recruits)));
                                NumberRecruitsPerSpecies(unique_species(i))=num_recruits;

                                 if num_recruits>0
                                    [xInd,yInd,zInd]=ind2sub(size(Recruits),find(Recruits>0));
                                    while num_recruits>size(xInd,1)
                                        Recruits(ind2sub(size(Recruits),find(Recruits>0)))=...
                                        Recruits(ind2sub(size(Recruits),find(Recruits>0)))-1;
                                        [xxInd,yyInd,zzInd]=ind2sub(size(Recruits),find(Recruits>0));
                                        xInd=[xInd;xxInd];
                                        yInd=[yInd;yyInd];
                                        zInd=[zInd;zzInd];
                                    end
                                    vec_recruits=length(E)+1:length(E)+size(xInd,1);

                                    %Copy species information to Epiphyte matrix
                                    E(vec_recruits,1:NumColSpeciesPool)= repmat(SpeciesPool(unique_species(i),:),size(xInd,1),1);
                                    E(vec_recruits,ColX)=xInd;
                                    E(vec_recruits,ColY)=yInd;
                                    E(vec_recruits,ColZ)=zInd;
                                    E(vec_recruits,ColMass)=0; %Initial size
                                    E(vec_recruits,ColStatus)=1; %status 1:alive
                                    E(vec_recruits,ColIndividualID)=(MaxIndividualID+1):((MaxIndividualID)+length(xInd)); %individual ID
                                    MaxIndividualID=MaxIndividualID+length(xInd);
                                 end
                            end
                       end
                    end

                    NumberRecruits=size(E(E(:,ColStatus)==1 ,:),1)-IntialNumberIndividualsTotal;
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                    E(E(:,1)==0,:)=[]; %in rare case, some individuals with only zeros are creates, which is wrong. This is to prevent the script to stop.
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % Growth
                    for i=1:size(E,1) %maybe it is faster if I do not use the if statement => speed testing
                        if(E(i,ColStatus)==1)
                            E(i,ColMass)=E(i,ColMass)+...
                            max(0,GrowthRate(E(i,ColMaximumMass),E(i,ColMass),E(i,ColGrowthRate))*...
                            Parabol(E(i,ColLightResponseA),E(i,ColLightResponseB),E(i,ColLightResponseC)...
                            ,Microhabitat(E(i,ColX),E(i,ColY),E(i,ColZ),3)));
                        end

                        %Add information about the voxel to the epiphyte matrix
                        E(i,ColSurfaceAreaOccupied)=(E(i,ColMass)^(2/3))/SurfaceBiomassScaling;
                        E(i,ColTotalSurfaceInVoxel)=Microhabitat(E(i,ColX),E(i,ColY),E(i,ColZ),1); %Total surface in voxel
                        E(i,ColSurfaceLossInVoxel)=Microhabitat(E(i,ColX),E(i,ColY),E(i,ColZ),2); %Percentage surface loss in this year
                        E(i,ColLightInVoxel)=Microhabitat(E(i,ColX),E(i,ColY),E(i,ColZ),3); %Light conditions in voxel

                    end
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %Mortality
                    for i=1:size(E,1)
                        if(E(i,ColStatus)==1)
                             %Mortality due to branch fall
                            if rand<Microhabitat(E(i,ColX),E(i,ColY),E(i,ColZ),2)
                                E(i,ColStatus)=3;

                            %Mortality due to changing light conditions
                            elseif  Microhabitat(E(i,ColX),E(i,ColY),E(i,ColZ),3)<E(i,ColMinLight) ||...
                                    Microhabitat(E(i,ColX),E(i,ColY),E(i,ColZ),3)>E(i,ColMaxLight)
                                E(i,ColStatus)=4;

                            %Natural mortality rate
                            elseif MortalityMethod==0 && rand<MortRateRandom
                                E(i,ColStatus)=5;

                            elseif MortalityMethod==1 && rand<(MortRateMass*(E(i,ColMass)^(MortRateMassScaling)))
                                E(i,ColStatus)=5;
                            end
                        end
                    end

                    %Mortality due to competition for space

                    %Calculate total surface area occupied by epiphytes per voxel
                    TotalSurfaceArePerVoxelOccupied=zeros(dimPlot(1),dimPlot(2),dimPlot(3));
                    for w=1:size(E,1)
                        if(E(w,ColStatus)==1)
                            TotalSurfaceArePerVoxelOccupied(E(w,ColX),E(w,ColY),E(w,ColZ))=TotalSurfaceArePerVoxelOccupied(E(w,ColX),E(w,ColY),E(w,ColZ))+E(w,ColSurfaceAreaOccupied);
                        end
                    end

                    %Indices of voxel where total area of epiphytes exeeds the available surface area
                    [IndX, IndY, IndZ]= ind2sub(size(TotalSurfaceArePerVoxelOccupied),find(TotalSurfaceArePerVoxelOccupied>Microhabitat(:,:,:,1)));

                    for i=1:size(IndX,1)

                        %Get all epis in voxel
                        EpisInVoxel=E(E(:,ColX)==IndX(i) & E(:,ColY)==IndY(i)...
                            & E(:,ColZ)==IndZ(i) & E(:,ColStatus)==1 ,:);

                        %Sort them by size (CompetitionMethod=1) or randomly (CompetitionMethod=2)
                        if CompetitionMethod==1
                             EpisInVoxel=sortrows(EpisInVoxel,-ColSurfaceAreaOccupied);
                        elseif CompetitionMethod==2
                             EpisInVoxel=EpisInVoxel(randperm(size(EpisInVoxel,1)),:);
                        end

                        CumulativeSumOfSurface=cumsum(EpisInVoxel(:,ColSurfaceAreaOccupied));
                        CumulativeSumOfSurfaceSum=sum(CumulativeSumOfSurface<=Microhabitat(IndX(i),IndY(i),IndZ(i),1));

                        if CumulativeSumOfSurfaceSum<size(EpisInVoxel,1)
                            E(ismember(E(:,ColIndividualID),EpisInVoxel(CumulativeSumOfSurfaceSum+1:size(EpisInVoxel,1),ColIndividualID)),ColStatus)=2;
                        end

                    end
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                    %Increase age
                    E(:,ColAge)=E(:,ColAge)+1;

                    %Save number of mortality event
                    MortalityCompetition=size(E(E(:,ColStatus)==2),1);
                    MortalityBranchFall=size(E(E(:,ColStatus)==3),1);
                    MortalityLight=size(E(E(:,ColStatus)==4),1);
                    MortalityNatural=size(E(E(:,ColStatus)==5),1);

                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %Store information in SummaryMatrixSpecies (summary over time for each species
                    for numSpecies=1:size(SpeciesPool,1)
                        SummaryMatrixSpecies(((numSpecies-1)*timeSteps)+t,ColSSpeciesID)=numSpecies;
                        SummaryMatrixSpecies(((numSpecies-1)*timeSteps)+t,ColSNumberIndividualsBeginning)=IntialNumberIndividuals(numSpecies);
                        SummaryMatrixSpecies(((numSpecies-1)*timeSteps)+t,ColSNumberIndividualsEnd)=sum(E(:,ColStatus)==1 & E(:,ColSpeciesID)==numSpecies);
                        SummaryMatrixSpecies(((numSpecies-1)*timeSteps)+t,ColSNumberMatureIndividuals)=sum(E(:,ColStatus)==1 & E(:,ColSpeciesID)==numSpecies & E(:,ColMass)>=E(:,ColMassAtMaturity));
                        SummaryMatrixSpecies(((numSpecies-1)*timeSteps)+t,ColSNumberRecruits)=NumberRecruitsPerSpecies(numSpecies);
                        SummaryMatrixSpecies(((numSpecies-1)*timeSteps)+t,ColSNumberMortalityBranchFall)=sum(E(:,ColStatus)==3 & E(:,ColSpeciesID)==numSpecies);
                        SummaryMatrixSpecies(((numSpecies-1)*timeSteps)+t,ColSNumberMortalityLight)=sum(E(:,ColStatus)==4 & E(:,ColSpeciesID)==numSpecies);
                        SummaryMatrixSpecies(((numSpecies-1)*timeSteps)+t,ColSNumberMortalityCompetition)=sum(E(:,ColStatus)==2 & E(:,ColSpeciesID)==numSpecies);
                        SummaryMatrixSpecies(((numSpecies-1)*timeSteps)+t,ColSNumberMortalityNatural)=sum(E(:,ColStatus)==5 & E(:,ColSpeciesID)==numSpecies);
                        if sum(E(:,ColStatus)==1 & E(:,ColSpeciesID)==numSpecies)&& IntialNumberIndividuals(numSpecies)>0
                            SummaryMatrixSpecies(((numSpecies-1)*timeSteps)+t,ColSNumberPopulationGrowthRate)=SummaryMatrixSpecies(((numSpecies-1)*timeSteps)+t,ColSNumberIndividualsEnd)/SummaryMatrixSpecies(((numSpecies-1)*timeSteps)+t,ColSNumberIndividualsBeginning);
                            SummaryMatrixSpecies(((numSpecies-1)*timeSteps)+t,ColSNumberPopulationGrowthRateLog)=log(SummaryMatrixSpecies(((numSpecies-1)*timeSteps)+t,ColSNumberPopulationGrowthRate));
                            SummaryMatrixSpecies(((numSpecies-1)*timeSteps)+t,ColSNumberBirthRate)=NumberRecruitsPerSpecies(numSpecies)/IntialNumberIndividuals(numSpecies);
                            SummaryMatrixSpecies(((numSpecies-1)*timeSteps)+t,ColSNumberDeathRate)=(sum(E(:,ColStatus)==3 & E(:,ColSpeciesID)==numSpecies)+sum(E(:,ColStatus)==4 & E(:,ColSpeciesID)==numSpecies)+sum(E(:,ColStatus)==2 & E(:,ColSpeciesID)==numSpecies)+sum(E(:,ColStatus)==5 & E(:,ColSpeciesID)==numSpecies))/IntialNumberIndividuals(numSpecies);
                            SummaryMatrixSpecies(((numSpecies-1)*timeSteps)+t,ColSAverageSize)=mean(E(E(:,ColSpeciesID)==numSpecies,ColMass));
                            SummaryMatrixSpecies(((numSpecies-1)*timeSteps)+t,ColSAverageAge)=mean(E(E(:,ColSpeciesID)==numSpecies,ColAge));
                            SummaryMatrixSpecies(((numSpecies-1)*timeSteps)+t,ColSMinLight)=min(E(E(:,ColSpeciesID)==numSpecies,ColLightInVoxel));
                            SummaryMatrixSpecies(((numSpecies-1)*timeSteps)+t,ColSMaxLight)=max(E(E(:,ColSpeciesID)==numSpecies,ColLightInVoxel));
                            SummaryMatrixSpecies(((numSpecies-1)*timeSteps)+t,ColSMeanLight)=mean(E(E(:,ColSpeciesID)==numSpecies,ColLightInVoxel));
                            SummaryMatrixSpecies(((numSpecies-1)*timeSteps)+t,ColSMinHeight)=min(E(E(:,ColSpeciesID)==numSpecies,ColZ));
                            SummaryMatrixSpecies(((numSpecies-1)*timeSteps)+t,ColSMaxHeight)=max(E(E(:,ColSpeciesID)==numSpecies,ColZ));
                            SummaryMatrixSpecies(((numSpecies-1)*timeSteps)+t,ColSMeanHeight)=mean(E(E(:,ColSpeciesID)==numSpecies,ColZ));
                        else
                            SummaryMatrixSpecies(((numSpecies-1)*timeSteps)+t,ColSNumberPopulationGrowthRate)=NaN;
                            SummaryMatrixSpecies(((numSpecies-1)*timeSteps)+t,ColSNumberPopulationGrowthRateLog)=NaN;
                            SummaryMatrixSpecies(((numSpecies-1)*timeSteps)+t,ColSNumberBirthRate)=NaN;
                            SummaryMatrixSpecies(((numSpecies-1)*timeSteps)+t,ColSNumberDeathRate)=NaN;
                            SummaryMatrixSpecies(((numSpecies-1)*timeSteps)+t,ColSAverageSize)=NaN;
                            SummaryMatrixSpecies(((numSpecies-1)*timeSteps)+t,ColSAverageAge)=NaN;
                            SummaryMatrixSpecies(((numSpecies-1)*timeSteps)+t,ColSMinLight)=NaN;
                            SummaryMatrixSpecies(((numSpecies-1)*timeSteps)+t,ColSMaxLight)=NaN;
                            SummaryMatrixSpecies(((numSpecies-1)*timeSteps)+t,ColSMeanLight)=NaN;
                            SummaryMatrixSpecies(((numSpecies-1)*timeSteps)+t,ColSMinHeight)=NaN;
                            SummaryMatrixSpecies(((numSpecies-1)*timeSteps)+t,ColSMaxHeight)=NaN;
                            SummaryMatrixSpecies(((numSpecies-1)*timeSteps)+t,ColSMeanHeight)=NaN;
                        end

                        SummaryMatrixSpeciesSave(((numSpecies-1)*timeSteps)+t,1)=t;
                    	SummaryMatrixSpeciesSave(((numSpecies-1)*timeSteps)+t,2:TotalColsSpeciesMatrix+1)=SummaryMatrixSpecies(((numSpecies-1)*timeSteps)+t,:);
                    end

                    %SummaryMatrixSpeciesSave((NumberOfSpecies*(t-1))+1:(NumberOfSpecies*t),1)=t;
                    %SummaryMatrixSpeciesSave((NumberOfSpecies*(t-1))+1:(NumberOfSpecies*t),2:TotalColsSpeciesMatrix+1)=SummaryMatrixSpecies;
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %Store information in SummaryMatrixCommunity
                    SummaryMatrixCommunity(t,1)=t; %TimeStep
                    SummaryMatrixCommunity(t,2)=InitialNumberSpecies; %NumberOfSpecies at beginning
                    SummaryMatrixCommunity(t,3)=length(unique(E(E(:,ColStatus)==1,ColSpeciesID))); %NumberOfSpecies at end
                    SummaryMatrixCommunity(t,4)=IntialNumberIndividualsTotal; %NumberIndividuals at beginning
                    SummaryMatrixCommunity(t,5)=size(E(E(:,ColStatus)==1 ,:),1); %NumberIndividuals at end
                    SummaryMatrixCommunity(t,6)=NumberRecruits; %Recruits
                    SummaryMatrixCommunity(t,7)=MortalityBranchFall; %MortalityBranchFall
                    SummaryMatrixCommunity(t,8)=MortalityLight; %MortalityLight
                    SummaryMatrixCommunity(t,9)=MortalityCompetition; %MortalityCompetition
                    SummaryMatrixCommunity(t,10)=MortalityNatural; %MortalityNatural
                    SummaryMatrixCommunity(t,11)=sum(sum(sum(Microhabitat(:,:,:,1))))/(dimPlot(1)*dimPlot(2)); %BranchSurfaceIndex
                    SummaryMatrixCommunity(t,12)=(sum(E(:,ColMass).^(2/3))/SurfaceBiomassScaling)/(sum(sum(sum(Microhabitat(:,:,:,1))))); %EpiphyteFilling
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                    %Command window information
                    fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% \n');
                    fprintf('Time step %d \n', t);
                    fprintf('Number of individuals %d \n', SummaryMatrixCommunity(t,5));
                    fprintf('Number of species %d \n', SummaryMatrixCommunity(t,3));
                    fprintf('Number of recruits %d \n', NumberRecruits);
                    fprintf('MortalityBranchFall %d \n', MortalityBranchFall);
                    fprintf('MortalityLight %d \n', MortalityLight);
                    fprintf('MortalityCompetition %d \n', MortalityCompetition);
                    fprintf('MortalityNatural %d \n', MortalityNatural);
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %Saving

                    %Save Epiphyte matrix for every time step
                    ESave=E(:,ColumsToSave);
                    EpiphyteMatrixCSV=strcat(DirectoryModelResultsRun,'\IndividualMatrixTimeStep',num2str(t),'.csv');
                    dlmwrite(EpiphyteMatrixCSV, ESave,'delimiter', '\t')

                    %Save SummaryMatrixSpecies for every time step
                    SummaryMatrixSpeciesCSV=strcat(DirectoryModelResultsRun,'\SpeciesSummary.csv');
                    dlmwrite(SummaryMatrixSpeciesCSV, SummaryMatrixSpeciesSave,'delimiter', '\t')

                    %Save SummaryMatrixCommunity for every time step (overwrite old one)
                    SummaryMatrixCommunityCSV=strcat(DirectoryModelResultsRun,'\CommunitySummary.csv');
                    dlmwrite(SummaryMatrixCommunityCSV, SummaryMatrixCommunity,'delimiter', '\t')
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                    %Remove dead individuals from Epimatrix
                    E(E(:,ColStatus)>1,:)=[];
                end
             end
        end
    end
end
