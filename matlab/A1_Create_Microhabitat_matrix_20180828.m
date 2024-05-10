%Create microhabitat matrices
clear all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Parameters that need to be specified/checked before running this script

%This parameter determines which type of microhatiat matrices are generated: 1: real GroIMP forest with dynamics
%2: static GroIMP forest (only forest at timeStepStart is used)
MicrohabitatType=1;

%Parameters of light model
kL=0.6; %light extinction coefficient
DistVoxToConsider=8; %How many ring around focal voxel to consider in light model (5 voxels in x and y direction)

%Choose the forest parameters that shall be calculated and stored in the microhabitat matrix (this list can be extended for possible new
%applications of the epiphyte model. 1: use this variable; 0: do not use it
TotalSurfaceAreaOpt=1;
SurfaceAreaLossOpt=1;
LightConditionsOpt=1;
AverageWeightedAngles=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Parameters that need to be specified if MicrohabitatType=1 or  MicrohabitatType=2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Directory of GroIMP files (this directory is stored in the Microhabitat folder so that the connection to the input GroIMP files is always clear)
DirectoryGroIMP='C:\Gunnar\ForestModel_20160129_BestModel_100x100_T1000';
%DirectorySaveMain='\\OMNISCIENTIA\Members\Gunnar\EpiphyteModels\MicrohabitatMatrices';
DirectorySaveMain='C:\Gunnar\EpiphyteModels\MicrohabitatMatrices';
DirectorySaveFolder='ForestModel_Best_100x100_T1000';

%Name under which the microhabitat matrices are saved (The name of the folder under
%which the microhabitat matrices are saved is standarized and only the name
%of the forest is required here, and only if MicrohabitatType=1 or MicrohabitatType=2)
NameForest='ForestModel_Best';
ReplicateForest=4;

%start and end timestep
timeStepStart=1; 
timeStepEnd=1000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Parameters that need to be specified if MicrohabitatType=3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%The following parameters are only needed if MicrohabitatType=3
%Dimensions of the theoretical forest
ForestHeight=40;
dimXTheoretical=50;
dimyTheoretical=50;
BAI=3; %branch area index for the static, theoretical forest
LAI=6; %leaf are index
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%The following parameters are generated automatically
%Load dimensions of forest patch from global forest file
if MicrohabitatType==1 || MicrohabitatType==2
    GlobalForest=dlmread(strcat(DirectoryGroIMP,'\Model\Forest_param_global.txt'), '\t',0,1);
    dimPlot=[GlobalForest(3),GlobalForest(4),GlobalForest(5)];
    corridor=GlobalForest(6);
    dimX=GlobalForest(3)+2*corridor; %MaxX+2*Corridor
    dimY=GlobalForest(4)+2*corridor; %MaxY+2*Corridor
    dimZ=GlobalForest(5); %MaxZ
end

%Create folder to save the microhabitat matrices
%The names of the folders are standadized:
%MicrohabitatType=1: 'Microhabitat_NameOfForestModel_SpatialExtent_timeSteps'
%MicrohabitatType=2: 'Microhabitat_NameOfForestModel_SpatialExtent_timeStep'
%MicrohabitatType=3: 'Microhabitat_BAI_LAI_kL'
if MicrohabitatType==1 
    NameMicrohabitatMatrix=strcat('Microhabitat_',NameForest,'_',num2str(dimPlot(1)),'x',num2str(dimPlot(2)),'x',num2str(dimPlot(3)),'_Rep',num2str(ReplicateForest));
    DirectoryMatrices=strcat(DirectorySaveMain,'\DynamicForests\',DirectorySaveFolder,'\',NameMicrohabitatMatrix);
elseif MicrohabitatType==2
    NameMicrohabitatMatrix=strcat('Microhabitat_',NameForest,'_',num2str(dimPlot(1)),'x',num2str(dimPlot(2)),'x',num2str(dimPlot(3)),'_Rep',num2str(ReplicateForest));
    DirectoryMatrices=strcat(DirectorySaveMain,'\StaticForests\',DirectorySaveFolder,'\',NameMicrohabitatMatrix);
elseif MicrohabitatType==3
    NameMicrohabitatMatrix=strcat('Microhabitat_BAI',num2str(BAI),'_LAI',num2str(LAI),'_kL',num2str(kL));
    DirectoryMatrices=strcat(DirectorySaveMain,'\UniformForests\',DirectorySaveFolder,'\',NameMicrohabitatMatrix);
end
mkdir(DirectoryMatrices)

%Copy global and pass forest file to microhabitat folder
if MicrohabitatType==1 || MicrohabitatType==2
    copyfile(strcat(DirectoryGroIMP,'\Model\Forest_param_global.txt'),DirectoryMatrices)
    copyfile(strcat(DirectoryGroIMP,'\Model\Forest_param_pass',num2str(ReplicateForest),'.txt'),DirectoryMatrices)
end

%Additional parameters
%Names of essential GroIMP files
shootFile=strcat('shoots_replicate_',num2str(ReplicateForest),'_time_step_');
trunkFile=strcat('trees_replicate_',num2str(ReplicateForest),'_time_step_');
voxelFile=strcat('voxel_replicate_',num2str(ReplicateForest),'_time_step_');

C = [0 0 1]; % Vector orthogonal to plane of X and Y
TotalVoxels=(DistVoxToConsider*2+1)^2;%Total number of adjacent voxels considered
MatrixDimension=sum([TotalSurfaceAreaOpt SurfaceAreaLossOpt LightConditionsOpt AverageWeightedAngles]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generation of microhabitat matrix of static or dynamic forest (MicrohabitatType=1 or MicrohabitatType=2)
% Here, the choosen parameter (total surface, surface loss, light conditions,average angle) are calculated for each voxel in each timestep

if MicrohabitatType==1 || MicrohabitatType==2
    
    %In a staic forest, only the initial forest at time step timeStepStart is of interest
    if MicrohabitatType==2
        timeStepEnd=timeStepStart+1;
    end
    
    for i=timeStepStart:timeStepEnd-1
        fprintf('Time step %d \n', i);
        tic

        %Load shoot and trunk files of actual and next timestep: Shoots at begin of year
        %and at the end of year/begin of next year
        if i~=timeStepStart
            ShootsBegin=ShootsEnd;
            ShootsEnd=dlmread(strcat(DirectoryGroIMP,'\Results\',shootFile,num2str(i+1),'.txt'),'\t',2,0);
            TrunksBegin=TrunksEnd;
            TrunksEnd=dlmread(strcat(DirectoryGroIMP,'\Results\',trunkFile,num2str(i+1),'.txt'),'\t',9,0);
        else
            ShootsBegin=dlmread(strcat(DirectoryGroIMP,'\Results\',shootFile,num2str(i),'.txt'),'\t',2,0);
            ShootsEnd=dlmread(strcat(DirectoryGroIMP,'\Results\',shootFile,num2str(i+1),'.txt'),'\t',2,0);
            TrunksBegin=dlmread(strcat(DirectoryGroIMP,'\Results\',trunkFile,num2str(i),'.txt'),'\t',9,0);
            TrunksEnd=dlmread(strcat(DirectoryGroIMP,'\Results\',trunkFile,num2str(i+1),'.txt'),'\t',9,0);
        end

        %inititalize all matrices
        Mat_surface_per_cell=zeros(dimX,dimY,dimZ);
        Mat_weighted_angle_per_cell=zeros(dimX,dimY,dimZ);
        Mat_light_per_cell=zeros(dimX,dimY,dimZ);
        Mat_surfaceloss_per_cell=zeros(dimX,dimY,dimZ);
        Mat_leafArea_per_cell=zeros(dimX,dimY,dimZ);
        
        %Get all branch segments that die during time step
        DeadSegments=ShootsBegin(~ismember(ShootsBegin(:,1),ShootsEnd(:,1)),1);
        [DeadSegments1, locDeadSegments]=ismember(DeadSegments,ShootsBegin(:,1));
        CounterDead=1;
        TotalDead=length(locDeadSegments);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Loop through all shoots and calculate total surface area, surface loss and weighted angles
        for j=1:length(ShootsBegin)

            %Get voxel the shoot is intersecting with
            %Here, we are using a simple approach because the shoots are
            %usually small only intersect with maximum  two voxels
            UniqueX=[ceil(ShootsBegin(j,8)) ceil(ShootsBegin(j,11))];
            UniqueY=[ceil(ShootsBegin(j,9)) ceil(ShootsBegin(j,12))];
            UniqueZ=[ceil(ShootsBegin(j,10)) ceil(ShootsBegin(j,13))];

            if UniqueX(1)==UniqueX(2)
                numX=1;
            else
                uniqueVectorX = UniqueX([true;diff(UniqueX(:))>0]);
                numX=length(uniqueVectorX);
            end

            if UniqueY(1)==UniqueY(2)
                numY=1;
            else
                uniqueVectorY = UniqueY([true;diff(UniqueY(:))>0]);
                numY=length(uniqueVectorY);
            end

            if UniqueZ(1)==UniqueZ(2)
                numZ=1;
            else
                uniqueVectorZ = UniqueZ([true;diff(UniqueZ(:))>0]);
                numZ=length(uniqueVectorZ);
            end

            for x=1:numX
                for y=1:numY
                    for z=1:numZ
                        
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        % Calculate new total surface area per voxel
                        if TotalSurfaceAreaOpt==1
                            Mat_surface_per_cell(UniqueX(x),UniqueY(y),UniqueZ(z))= ...
                                Mat_surface_per_cell(UniqueX(x),UniqueY(y),UniqueZ(z))+...
                                ((ShootsBegin(j,5)/(numX*numY*numZ))*ShootsBegin(j,6)*pi/2); %Surface of single branch within voxel
                        end
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        % If branch is lost during this time step, add it
                        % to lost surface
                        if SurfaceAreaLossOpt==1
                            if(j==locDeadSegments(CounterDead))
                                %fprintf('Dead branch segment %d \n',CounterDead);
                                %test=test+1;
                                CounterDead=min(TotalDead,CounterDead+1);

                                Mat_surfaceloss_per_cell(UniqueX(x),UniqueY(y),UniqueZ(z))= ...
                                    Mat_surfaceloss_per_cell(UniqueX(x),UniqueY(y),UniqueZ(z))+...
                                    ((ShootsBegin(j,5)/(numX*numY*numZ))*ShootsBegin(j,6)*pi/2);

                            end
                        end
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        % Calculate weighted angles per voxel
                        if AverageWeightedAngles==1
                            % Calculate angle for the shoot relative to the plane x
                            % and y
                            position1=[ShootsBegin(j,8), ShootsBegin(j,9),ShootsBegin(j,10)];
                            position2=[ShootsBegin(j,11), ShootsBegin(j,12),ShootsBegin(j,13)];
                            V=position2-position1;
                            alpha=sum(C.*V)/(sqrt(V(1)^2+V(2)^2+V(3)^2)*sqrt(C(1)^2+C(2)^2+C(3)^2));
                            ShootAngle=abs(90-(acos(alpha)/pi*180)); %Angle for shoot

                            % Calculate weighted angle for the voxel
                            Mat_weighted_angle_per_cell(UniqueX(x),UniqueY(y),UniqueZ(z))= ...
                                (...
                                ((Mat_surface_per_cell(UniqueX(x),UniqueY(y),UniqueZ(z)))-((ShootsBegin(j,5)/(numX*numY*numZ))*ShootsBegin(j,6)*pi/2))/...
                                Mat_surface_per_cell(UniqueX(x),UniqueY(y),UniqueZ(z))*...
                                Mat_weighted_angle_per_cell(UniqueX(x),UniqueY(y),UniqueZ(z)))...
                                +(...
                                 (((ShootsBegin(j,5)/(numX*numY*numZ))*ShootsBegin(j,6)*pi/2))/...
                                 Mat_surface_per_cell(UniqueX(x),UniqueY(y),UniqueZ(z))*...
                                ShootAngle);
                        end
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        
                    end
                end
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Loop through all trunks and calculate total surface area, surface loss and weighted angles
        
        %Get all trees that die during time step
        DeadSegments=TrunksBegin(~ismember(TrunksBegin(:,1),TrunksEnd(:,1)),1);
        [DeadSegments1, locDeadSegments]=ismember(DeadSegments,TrunksBegin(:,1));
        CounterDead=1;
        TotalDead=length(locDeadSegments);

        for j=1:length(TrunksBegin)

                X=ceil(TrunksBegin(j,6)); %X voxel of tree
                Y=ceil(TrunksBegin(j,7)); %Y voxel of tree
                Height=TrunksBegin(j,3); %Height of tree
                Diameter=TrunksBegin(j,4); %Diameter of tree

                SurfaceAreaTotal=0;

            for Z=ceil(Height):-1:1

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Calculate new total surface area per voxel (for trunks)
                if TotalSurfaceAreaOpt==1
                    hCone=Height-Z+1; %height of cylinder from top to bottom of voxel
                    %rCone=(hCone/Height)*(Diameter/2); %radius of cylinder at bottom of voxel
                    rCone=(Diameter/2); %radius of cylinder at bottom of voxel

                    % Calculate total surface area in voxel and save it in matrix
                    SurfaceAreaInVoxel=pi*rCone*sqrt((rCone^2)+(hCone^2))-SurfaceAreaTotal;
                    Mat_surface_per_cell(X,Y,Z)=Mat_surface_per_cell(X,Y,Z)+SurfaceAreaInVoxel;

                    % Update total surface area of cylinder so far (to use in next step)
                    SurfaceAreaTotal=SurfaceAreaInVoxel+SurfaceAreaTotal;
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % If trunk is lost during this time step, add it
                % to lost surface
                if SurfaceAreaLossOpt==1
                    if(j==locDeadSegments(CounterDead))
                        CounterDead=min(TotalDead,CounterDead+1);
                        Mat_surfaceloss_per_cell(X,Y,Z)=Mat_surfaceloss_per_cell(X,Y,Z)+...
                            SurfaceAreaInVoxel;
                    end
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Update weighted angle for the voxel
                if AverageWeightedAngles==1
                    Mat_weighted_angle_per_cell(X,Y,Z)= ...
                        ( (Mat_surface_per_cell(X,Y,Z)-SurfaceAreaInVoxel)/...
                        Mat_surface_per_cell(X,Y,Z)*Mat_weighted_angle_per_cell(X,Y,Z))...
                        +((SurfaceAreaInVoxel/Mat_surface_per_cell(X,Y,Z)*90)); % upright 90° angle assumed
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%         

            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
        %Calculate light conditions in voxels (relative light conditions)
        if LightConditionsOpt==1
        
            %Load file containing information about leaf area per voxel
            Voxels=dlmread(strcat(DirectoryGroIMP,'\Results\',voxelFile,num2str(i),'.txt'),'\t',2,0);
            Voxels(:,1:3)=Voxels(:,1:3)+1; %Voxel file start with x=y=z=0 => synchronize with matrices used here

            %Store information on leaf area in matrix
            for j=1:length(Voxels)
                  Mat_leafArea_per_cell(Voxels(j,1),Voxels(j,2),Voxels(j,3))=Voxels(j,4);
            end

            %Calculate single column light conditions based on leaf area distribution
            for x=1:dimX
                for y=1:dimY
                    for z=1:dimZ

                         Mat_light_per_cell(x,y,z)=exp(-kL*(sum(Mat_leafArea_per_cell(x,y,z:dimZ))/10000));

                    end
                end
            end

            %Copy light conditions
            Mat_light_per_cell_Copy=Mat_light_per_cell;

            %Calculate final light conditions by accounting for the light
            %conditions in adjacent voxels
            for x=corridor:dimX-corridor
                 for y=corridor:dimY-corridor
                    for z=1:dimZ

                        TotalContribution=0;

                        %loop over ring surrounding the focal voxel
                        for xx=x-DistVoxToConsider:x+DistVoxToConsider
                           for yy=y-DistVoxToConsider:y+DistVoxToConsider

                               Ring=max(abs(xx-x),abs(yy-y));
                               Contribution=(1/(DistVoxToConsider+1))*(1/max(1,(Ring*8)))*Mat_light_per_cell_Copy(xx,yy,z);
                               TotalContribution=TotalContribution+Contribution;

                           end
                       end

                        Mat_light_per_cell(x,y,z)=TotalContribution;

                    end
                end
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
        %Store information in Microhabitat matrix and save matrix for this
        %timestep
        Microhabitat=zeros(dimPlot(1),dimPlot(2),dimPlot(3),MatrixDimension,'single');%eventuell nur 5 dimensionen um platz zu sparen
        
        if TotalSurfaceAreaOpt==1
            Microhabitat(:,:,:,1)=Mat_surface_per_cell(corridor+1:dimX-corridor,corridor+1:dimY-corridor,1:dimZ);
        end
        
        if SurfaceAreaLossOpt==1
            if MicrohabitatType==1
                Microhabitat(:,:,:,2)=Mat_surfaceloss_per_cell(corridor+1:dimX-corridor,corridor+1:dimY-corridor,1:dimZ)./Mat_surface_per_cell(corridor+1:dimX-corridor,corridor+1:dimY-corridor,1:dimZ);
            end
            
            if MicrohabitatType==2
                Microhabitat(:,:,:,2)=0;
            end
        end
        
        if LightConditionsOpt==1
            Microhabitat(:,:,:,3)=Mat_light_per_cell(corridor+1:dimX-corridor,corridor+1:dimY-corridor,1:dimZ);
        end
        
        if AverageWeightedAngles==1
            Microhabitat(:,:,:,4)=Mat_weighted_angle_per_cell(corridor+1:dimX-corridor,corridor+1:dimY-corridor,1:dimZ);
        end
               
        MicrohabitatMatSave=strcat(DirectoryMatrices,'\MicrohabitatMatrix',num2str(i),'.mat');
        save(MicrohabitatMatSave,'Microhabitat')
        toc 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    end
    
    %Save the directory of the GroIMP file
    dlmwrite(strcat(DirectoryMatrices,'\DirectoryForestModelGroIMP.txt'),DirectoryGroIMP,'') %Save associated directory of species pools

    %Save dimensions of plot in seperate file
    save(strcat(DirectoryMatrices,'\dimPlot.mat'),'dimPlot')
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
