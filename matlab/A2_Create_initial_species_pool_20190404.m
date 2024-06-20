%Create species matrices
clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Parameters that need to be specified/checked before running this script

%Folder to save species trait matrices
NameSpeciesPool='IntAgeMat_2_IntRec_70'; %Give meaningful name (the species type is automatically added to the name)

%Define number of species in species pool and total number of species pools to be created
numSpeciesPools=100;
NumberOfSpecies=100;

%The following option defines if correlations between traits are consider or not
CorrelationMassAgeOfMaturity=1; %Correlation between the mass and the age of maturity (this also influences the growth rate)
CorrelationMassRecruitment=1; %Correlation between the mass and the recruitment

InterceptAgeMaturity=2;
ScalingAgeMaturity=0.25; %Scaling factor according to metabolic theory
AgeMaturityMetabolic=@(InterceptAgeMaturity,ScalingAgeMaturity,Mass) (InterceptAgeMaturity*(Mass^ScalingAgeMaturity));


%If correlations are choosen, the following parameters define the shape of the correlations
%1. Correlations if CorrelationMassAgeOfMaturity=1
MaxMassRangeCorr=[2 3000]; %maximum mass of species/functional types (g)
%%%AgeAtMaturityRangeCorr=[1 15]; %age at which maturity is reaches (years) #Comment2019 => should not ne needed, delete
AgeAtMaturityDevCorr=0.25; %relative deviation from mean age of maturity

%2. Correlations if CorrelationMassRecruitment=1
RecruitmentNormalizeAtSize1Corr=70; %Factor converting the reproductive biomass to potential recruits
SlopeRecruitmentCorr=0*RecruitmentNormalizeAtSize1Corr; %Slope of the correlation between mass and recruitment
RecruitmentInvestmentRelMeanCorr=0.1; %Anual investment in reproduction in relation to vegetative biomass (decrease due to correlation with mass)
RecruitmentInvestmentRelDevCorr=0.25; %The relative deviation from the mean recruitment
RecruitmentIncMaxCorr=0;

%Parameters of light model (needed to convert the height nicht and the light niche, these values do not have to
%be the same as used in the microhabitat matrices)
kL=0.6; %light extinction coefficient
Imax=900; %maximum light intensity
LAI=6; %leaf area index

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Define trait (ranges) if random species pool(SpeciesPoolType=0) is choose
%If no correlations between traits are choosen (CorrelationMassAgeOfMaturity=0 || CorrelationMassRecruitment=0), traits are randomly choosen from the following ranges
MaxMassRandom=[2 3000]; %maximum mass of species/functional types (g)
MaxMassLogScaleRandom=1; %define if the mass is choosen based on the log scale (MaxMassLogScale=1) or on the normal scale (MaxMassLogScale=0)
AgeAtMaturityRandom=[1 1]; %age at which maturity is reaches (years)
RecruitmentNormalizeAtSize1Random=[1 20]; %This parameter regulates the range of recruitment in thise cases
RecruitmentInvestmentRelMeanRandom=[0.07 0.12]; %Not that the effective recruitment is RecruitmentNormalizeAtSize1Random*RecruitmentInvestmentRelMeanRandom
RecruitmentIncRandom=[0 0];
MassAtMaturityRelativeRandom=[0.5 0.7]; %Relative mass in relation to maximum Size
HeightBreadthRandom=[0.15 0.7]; %Relative height
DispersalKernelRandom=[0.03 0.5]; %The higher this values, the more local is the dispersal
DispersalKernelAsymmetryRandom=[0.5 0.95]; %The trait describes the relative proportion of seed dispersed below the mother (i.e. 0.5=> symmetric dispersal kernel)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Define which trait is varied if a sequential species pool(SpeciesPoolType=1) is choose.
%For proper results, only one of the following traits should be defined as sequence, while for the other traits, invariable trait values should to be specified.
%If no correlations between traits are choosen (CorrelationMassAgeOfMaturity=0 || CorrelationMassRecruitment=0),the following traits are used
MaxMassSeq=[2 3000]; %maximum mass of species/functional types (g)
MaxMassLogScaleSeq=1; %define if the mass is choosen based on the log scale (MaxMassLogScale=1) or on the normal scale (MaxMassLogScale=0)
AgeAtMaturitySeq=2; %if no correlation is defined, this value is used
RecruitmentNormalizeAtSize1Seq=20; %This parameter regulates the range of recruitment in thise cases
RecruitmentInvestmentRelMeanSeq=0.1; %anual investment in reproduction in relation to vegetative biomass
RecruitmentIncSeq=0; %Increase in realtive reproductive allocation with mass 0: no increase; 1: doubling
MassAtMaturityRelativeSeq=0.5; %Relative mass in relation to maximum Size
HeightBreadthSeq=1; %Relative height
MeanHeightSeq=0.5;
DispersalKernelSeq=0; %The higher this values, the more local is the dispersal
DispersalKernelAsymmetrySeq=0.5; %The trait describes the relative proportion of seed dispersed below the mother (i.e. 0.5=> symmetric dispersal kernel)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Define traits if neutral species pool(SpeciesPoolType=2) is choose
%The traits are the same for all species and have to be specified below. If correlations between traits are choosen,
%the age at maturity and the recruitment are based on the correlations defined above instead of the one defined below
MaxMassNeutral=100; %maximum mass of species/functional types (g)
AgeAtMaturityNeutral=3; %age at which maturity is reaches (years)
RecruitmentNormalizeAtSize1Neutral=15; %If no correlation is choose, this value is used as recruitment (it is not multiplied by RecruitmentInvestmentRel!)
RecruitmentInvestmentRelMeanNeutral=0.1;
RecruitmentIncNeutral=0; %Increase in realtive reproductive allocation with mass 0: no increase; 1: doubling
MassAtMaturityRelativeNeutral=0.5; %Relative mass in relation to maximum Size
HeightBreadthNeutral=0.5; %Relative height
MeanHeightNeutral=0.5;  %Mean height
DispersalKernelNeutral=1.5; %The higher this values, the more local is the dispersal
DispersalKernelAsymmetryNeutral=0.5; %The trait describes the relative proportion of seed dispersed below the mother (i.e. 0.5=> symmetric dispersal kernel)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Create folder to save species trait matrices
if SpeciesPoolType==0
    if CorrelationMassAgeOfMaturity==1 || CorrelationMassRecruitment==1
        FullNameSpeciesPool=strcat('SP_Random_',NameSpeciesPool,'_TraitCorrOn');
    else
        FullNameSpeciesPool=strcat('SP_Random_',NameSpeciesPool,'_TraitCorrOff');
    end
elseif SpeciesPoolType==1
    if CorrelationMassAgeOfMaturity==1 || CorrelationMassRecruitment==1
        FullNameSpeciesPool=strcat('SP_Sequential_',NameSpeciesPool,'_TraitCorrOn');
    else
        FullNameSpeciesPool=strcat('SP_Sequential_',NameSpeciesPool,'_TraitCorrOff');
    end
elseif SpeciesPoolType==2
    if CorrelationMassAgeOfMaturity==1 || CorrelationMassRecruitment==1
        FullNameSpeciesPool=strcat('SP_Neutral_',NameSpeciesPool,'_TraitCorrOn');
    else
        FullNameSpeciesPool=strcat('SP_Neutral_',NameSpeciesPool,'_TraitCorrOff');
    end
end

SaveDirectory=strcat('\\OMNISCIENTIA\Members\Gunnar\EpiphyteModels\SpeciesPools','\',FullNameSpeciesPool);
mkdir(SaveDirectory);

%Copy this matlab file to the species pool folder (to see how the species
%pool was generated)
FileName=strcat(mfilename,'.m');
FileNameAndLocation=mfilename('fullpath');
currentfile=strcat(FileNameAndLocation, '.m');
newfile=strcat(SaveDirectory,'\',FileName);
copyfile(currentfile,newfile);

ColumnHeaders={'SpeciesID','MaximumMass','MassAtMaturity','GrowthRate',...
        'DispersalKernel','DispersalKernelAsymmetry','RecruitmentInvestmentRel','RecruitmentInc','MinLight','MaxLight','OptimumLight','LightBreadth',...
        'LightResponseA','LightResponseB','LightResponseC','MinHeightRel','MaxHeightRel',...
        'MeanHeightRel','HeightBreadth','MaxRecruitsAtMaxMass','MaxRecruitsAtMassAtMaturity','AgeAtMaturity'};

%Save ColumnHeader for further use
xlswrite(strcat(SaveDirectory,'\ColumnHeaders.xls'),ColumnHeaders);

%Save SpeciesPoolType for further use
save(strcat(SaveDirectory,'\SpeciesPoolType.mat'),'SpeciesPoolType')

%Save number of species pools and number of species for further use
save(strcat(SaveDirectory,'\numSpeciesPools.mat'),'numSpeciesPools')
save(strcat(SaveDirectory,'\NumberOfSpecies.mat'),'NumberOfSpecies')

%Save type of correlations
save(strcat(SaveDirectory,'\CorrelationMassAgeOfMaturity.mat'),'CorrelationMassAgeOfMaturity')
save(strcat(SaveDirectory,'\CorrelationMassRecruitment.mat'),'CorrelationMassRecruitment')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Main loop (for random generation of species pool)
for Num=1:numSpeciesPools

    %Trait matrix where the trait information of each species is saved
    SpeciesTraitMatrix=zeros(NumberOfSpecies,16);

    for NumSpecies=1:NumberOfSpecies

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Maximum Size, size at maturity and growth rate (here, we are
        % choosing from log scale because usually, there are more smaller species than larger ones)
        if MaxMassLogScaleRandom==1
            MaxMassLog=random('unif',log10(MaxMassRandom(1)),log10(MaxMassRandom(2)));
            MaxMass=10^MaxMassLog;
        elseif MaxMassLogScaleRandom==0
            MaxMass=random('unif',MaxMassRandom(1),MaxMassRandom(2));
        end

        %The age at maturity is the given by
        AgeAtMaturity=AgeMaturityMetabolic(InterceptAgeMaturity,ScalingAgeMaturity,MaxMass);
        AgeAtMaturity=AgeAtMaturity*random('unif',1-AgeAtMaturityDevCorr,1+AgeAtMaturityDevCorr); %Add stochasticity:

        %We assume that the mass at maturity is a function of the maximum size
        MassAtMaturity=random('unif',MassAtMaturityRelativeRandom(1),MassAtMaturityRelativeRandom(2))*MaxMass;

        %In the model, we are approximating growth by a Betalanffy growth curve, which generall is as follows
        %SizeFunctionOfAge=@(MaxMass,K,Age) (MaxMass*(1-exp(-K*(Age))));
        %By assuming that the Betalanffy growth curve crosses the point AgeAtMaturity/Size MaturityMassAtMaturity,
        %the growth rate K of this function can be calculated:
        K=-(log(1)+log(1-(MassAtMaturity/MaxMass)))/AgeAtMaturity;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Recruitment traits

        %2. Correlations if CorrelationMassRecruitment=1
        RecruitmentNormalizeAtSize1Corr=70; %Factor converting the reproductive biomass to potential recruits
        SlopeRecruitmentCorr=0*RecruitmentNormalizeAtSize1Corr; %Slope of the correlation between mass and recruitment
        RecruitmentInvestmentRelMeanCorr=0.1; %Anual investment in reproduction in relation to vegetative biomass (decrease due to correlation with mass)
        RecruitmentInvestmentRelDevCorr=0.25; %The relative deviation from the mean recruitment
        RecruitmentIncMaxCorr=0;

        if CorrelationMassRecruitment==1
            RecruitmentInvestmentRel=random('unif',RecruitmentInvestmentRelMeanCorr*(1-RecruitmentInvestmentRelDevCorr),RecruitmentInvestmentRelMeanCorr*(1+RecruitmentInvestmentRelDevCorr));
            RecruitmentNormalizeAtSize1=RecruitmentNormalizeAtSize1Corr; %Factor converting the reproductive biomass to potential recruits
            SlopeRecruitment=0; %Slope of the correlation between mass and recruitment
            InterceptRecruitment=RecruitmentNormalizeAtSize1; %
            RecruitmentInc=0;
        elseif CorrelationMassRecruitment==0
            RecruitmentInvestmentRel=random('unif',RecruitmentInvestmentRelMeanRandom(1),RecruitmentInvestmentRelMeanRandom(2));
            RecruitmentNormalizeAtSize1=random('unif',RecruitmentNormalizeAtSize1Random(1),RecruitmentNormalizeAtSize1Random(2));
            SlopeRecruitment=0; %No slope if no correlation is choosen
            InterceptRecruitment=RecruitmentNormalizeAtSize1-SlopeRecruitment;
            RecruitmentInc=random('unif',RecruitmentIncRandom(1),RecruitmentIncRandom(2)); %Not meaningful if no correlation
         end

        DispersalKernel=random('unif',DispersalKernelRandom(1),DispersalKernelRandom(2));
        DispersalKernelAsymmetry=random('unif',DispersalKernelAsymmetryRandom(1),DispersalKernelAsymmetryRandom(2));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Traits of ecologcial niche
        %1. Randomly choose mean height and height breadth
        MeanHeight=random('unif',0,1); %realtive height in relation to canopy height
        HeightBreadthTheoretical=random('unif',HeightBreadthRandom(1),HeightBreadthRandom(2));

        %Minimum and maximum height under which the species is
        %able to survive
        MinHeight=max(0,MeanHeight-(HeightBreadthTheoretical/2));
        MaxHeight=min(1,MeanHeight+(HeightBreadthTheoretical/2));
        HeightBreadth=MaxHeight-MinHeight;

        %2. Convert heigth ranges to light ranges
        %For this, a standard forest with the following parameters is
        %assumed
        MinLight=Imax*exp(-kL*LAI*(1-MinHeight));
        MaxLight=Imax*exp(-kL*LAI*(1-MaxHeight));
        OptimumLight=(MaxLight+MinLight)/2;
        LightBreadth=MaxLight-MinLight;

        %3. Calculate parameters of parabolic response curve y=ax^2+bx+c
        %We assume that the function is a paraboloid which goes trough
        %three points (MinLight/0) (MaxLight/0) (OptimumLight/1)

        x1=MinLight;
        x2=MaxLight;
        x3=OptimumLight;

        y1=0;
        y2=0;
        y3=1;

        a = (x1*(y2-y3)+x2*(y3-y1)+x3*(y1-y2))/((x1-x2)*(x1-x3)*(x3-x2));
        b = (x1^2*(y2-y3)+x2^2*(y3-y1)+x3^2*(y1-y2))/((x1-x2)*(x1-x3)*(x2-x3));
        c = (x1^2*(x2*y3-x3*y2)+x1*(x3^2*y2-x2^2*y3)+x2*x3*y1*(x2-x3))/((x1-x2)*(x1-x3)*(x2-x3));
        %In the model, based on these parameters the light response for each species can be calculated:
        %Parabol=@(a,b,c,x) a*x^2+b*x+c;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Assign trait values for each species
        SpeciesTraitMatrix(NumSpecies,1)=NumSpecies; %Species or functional type
        SpeciesTraitMatrix(NumSpecies,2)=MaxMass; %Maximum mass
        SpeciesTraitMatrix(NumSpecies,3)=MassAtMaturity; %Mass at maturity
        SpeciesTraitMatrix(NumSpecies,4)=K; %Species-specfic growth rate
        SpeciesTraitMatrix(NumSpecies,5)=DispersalKernel; %Dispersal: factor b in negative exp funtion
        SpeciesTraitMatrix(NumSpecies,6)=DispersalKernelAsymmetry; %Dispersal: factor b in negative exp funtion
        SpeciesTraitMatrix(NumSpecies,7)=RecruitmentInvestmentRel; %anual reproductive allocation in relation to vegetative biomass
        SpeciesTraitMatrix(NumSpecies,8)=RecruitmentInc; %Increase in realtive reproductive allocation with mass 0: no increase; 1: doubling
        SpeciesTraitMatrix(NumSpecies,9)=MinLight; %Min Light conditions
        SpeciesTraitMatrix(NumSpecies,10)=MaxLight; %Max Light conditions
        SpeciesTraitMatrix(NumSpecies,11)=OptimumLight; %Optimum Light conditions
        SpeciesTraitMatrix(NumSpecies,12)=LightBreadth; %Realised Light breadth
        SpeciesTraitMatrix(NumSpecies,13)=a; %Factor a of light response function
        SpeciesTraitMatrix(NumSpecies,14)=b; %Factor b of light response function
        SpeciesTraitMatrix(NumSpecies,15)=c; %Factor c of light response function
        SpeciesTraitMatrix(NumSpecies,16)=MinHeight; %Relative minimum height in a uniform standard forest
        SpeciesTraitMatrix(NumSpecies,17)=MaxHeight; %Relative maximum height in a uniform standard forest
        SpeciesTraitMatrix(NumSpecies,18)=MeanHeight; %Relative optimum height in a uniform standard forest
        SpeciesTraitMatrix(NumSpecies,19)=HeightBreadth; %Niche Breadth Height
        SpeciesTraitMatrix(NumSpecies,20)=(InterceptRecruitment)*RecruitmentInvestmentRel; %Potential maximum number of recruits at maximum mass
        SpeciesTraitMatrix(NumSpecies,21)=(InterceptRecruitment+SlopeRecruitment*MassAtMaturity)*RecruitmentInvestmentRel; %Potential maximum number of recruits at mass at maturity
        SpeciesTraitMatrix(NumSpecies,22)=AgeAtMaturity; %Average age at maturity under optimal conditions
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Save trait matrices
    %Create dataset from table (including headers)
    SpeciesTable=dataset({SpeciesTraitMatrix,ColumnHeaders{:}});

    %Save trait matrix
    SaveFile=strcat(SaveDirectory,'\SpeciesPool',num2str(Num),'.csv');
    dlmwrite(SaveFile,SpeciesTraitMatrix,'\t');

    SaveFileDataset=strcat(SaveDirectory,'\SpeciesPoolHeader',num2str(Num),'.txt');
    export(SpeciesTable,'file',SaveFileDataset)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Save trait ranges used to generate the species pool
    TraitRanges=zeros(20,2);
    TraitRanges(1,:)=SlopeRecruitment;
    TraitRanges(2,:)=InterceptRecruitment;
    TraitRanges(3,:)=MaxMassRangeCorr(1:length(MaxMassRangeCorr));
    %%TraitRanges(4,:)=AgeAtMaturityRangeCorr(1:length(AgeAtMaturityRangeCorr));  #Comment2019 => should not ne needed, delete
    TraitRanges(5,:)=AgeAtMaturityDevCorr(1:length(AgeAtMaturityDevCorr));
    TraitRanges(6,:)=RecruitmentNormalizeAtSize1Corr(1:length(RecruitmentNormalizeAtSize1Corr));
    TraitRanges(7,:)=SlopeRecruitmentCorr(1:length(SlopeRecruitmentCorr));
    TraitRanges(8,:)=RecruitmentInvestmentRelMeanCorr(1:length(RecruitmentInvestmentRelMeanCorr));
    TraitRanges(9,:)=RecruitmentInvestmentRelDevCorr(1:length(RecruitmentInvestmentRelDevCorr));
    TraitRanges(10,:)=RecruitmentIncMaxCorr(1:length(RecruitmentIncMaxCorr));
    TraitRanges(11,:)=MaxMassRandom(1:length(MaxMassRandom));
    TraitRanges(12,:)=AgeAtMaturityRandom(1:length(AgeAtMaturityRandom));
    TraitRanges(13,:)=RecruitmentNormalizeAtSize1Random(1:length(RecruitmentNormalizeAtSize1Random));
    TraitRanges(14,:)=RecruitmentInvestmentRelMeanRandom(1:length(RecruitmentInvestmentRelMeanRandom));
    TraitRanges(15,:)=RecruitmentIncRandom(1:length(RecruitmentIncRandom));
    TraitRanges(16,:)=MassAtMaturityRelativeRandom(1:length(MassAtMaturityRelativeRandom));
    TraitRanges(17,:)=HeightBreadthRandom(1:length(HeightBreadthRandom));
    TraitRanges(18,:)=DispersalKernelRandom(1:length(DispersalKernelRandom));
    TraitRanges(19,:)=DispersalKernelAsymmetryRandom(1:length(DispersalKernelAsymmetryRandom));
    TraitRanges(20,:)=MaxMassLogScaleRandom(1:length(MaxMassLogScaleRandom));
    dlmwrite(strcat(SaveDirectory,'\TraitRanges.csv'),TraitRanges,'\t');
end
