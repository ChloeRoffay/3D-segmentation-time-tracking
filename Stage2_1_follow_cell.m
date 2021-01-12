
function [MatrixCentroid_final_t,zMilieu_t] = Stage2_1_follow_cell(pathMovie,nameMovie,tmin,tmax,scale1D,nombreCelluleEtudie,zMilieu,k_animal,zmax)

%% Initialization

nombreCelluleMaxAttendues = 1000;
MatrixCentroid_brut = NaN(nombreCelluleMaxAttendues,2,tmax);
MatrixCentroid_t_brut = NaN(nombreCelluleMaxAttendues,2,tmax);

% Create dossier Data
status = mkdir([pathMovie filesep 'Data']);

%% Identify zMilieu

zMilieu_t = ones(tmax,1)*zMilieu;

%% Extract the centroid of all the objects for every t

for t=tmin:tmax
    
    BW = imread([pathMovie filesep 't' num2str(t,'%04d') filesep 'Output_results' filesep 'Unionseg_' nameMovie '_t' num2str(t,'%04d') '_z' num2str(zMilieu_t(t,1),'%04d') '.png']);
    
    image_CC = bwconncomp(BW,4);
    Stat_All_Regions = regionprops(image_CC,'Centroid');
    
    Stat_All_Regions_cell = (struct2cell(Stat_All_Regions))';
    %matCentroid = cell2mat(Stat_All_Regions_cell(:,1))*scale1D; % VRAIE ECHELLE
    matCentroid = cell2mat(Stat_All_Regions_cell(:,1)); % % ECHELLE MATLAB
    
    MatrixCentroid_brut(1:size(matCentroid,1),1,t) = matCentroid(:,1);
    MatrixCentroid_brut(1:size(matCentroid,1),2,t) = matCentroid(:,2);
    
end


%% Extract distance of all the centroids compared to tmin 
% Possibility to compared between t-1 and t but in my case the complexity was useless

matrice_R_brut = NaN(nombreCelluleMaxAttendues,nombreCelluleMaxAttendues,tmax);

for iteration_temps = tmin+1:tmax
    for iteration_cellule = 1:1:size(MatrixCentroid_brut,1)
        
        xA2 = MatrixCentroid_brut(:,1,iteration_temps);
        yA2 = MatrixCentroid_brut(:,2,iteration_temps);
        
        temps_initial = 1;
        xA1 = MatrixCentroid_brut(iteration_cellule,1,temps_initial);
        yA1 = MatrixCentroid_brut(iteration_cellule,2,temps_initial);
        
        r_carre = power((xA2-xA1),2)+power((yA2-yA1),2);
        r = power(r_carre,0.5);
        
        matrice_R_brut(1:size(r,1),iteration_cellule,iteration_temps) = r;
        
    end
end

%% Determine the closest centroid for each object of the image

% Select the minimum
min_matrice_R_brut = min(matrice_R_brut);

% Delete NaN
min_matrice_R_brut_tf = ~isnan(min_matrice_R_brut);

% Line number gives which centroid are connected (column number)
min_matrice_R_brut_tri = min_matrice_R_brut(min_matrice_R_brut_tf);
tf_brut = ismember(matrice_R_brut,min_matrice_R_brut_tri); 

%% Matrix is ordered accordingly

for iteration_temps = tmin+1:tmax
    tf_temps = tf_brut(:,:,iteration_temps);
    [row,col,~] = find(tf_temps>0);
    
    for iteration_cellule = 1:1:size(col,1)
        % x car seconde dimension = 1
        MatrixCentroid_t_brut(col(iteration_cellule),1,iteration_temps)= MatrixCentroid_brut(row(iteration_cellule),1,iteration_temps);
        
        % y car seconde dimension = 2
        MatrixCentroid_t_brut(col(iteration_cellule),2,iteration_temps)= MatrixCentroid_brut(row(iteration_cellule),2,iteration_temps);
    end
    
end

% Integrate t1 (since it is the comparaison it has to be added)
MatrixCentroid_t_brut(:,1,1)= MatrixCentroid_brut(:,1,1);
MatrixCentroid_t_brut(:,2,1)= MatrixCentroid_brut(:,2,1);

%% OSEF = delete the NaN

TF_nan = isnan(MatrixCentroid_t_brut(:,2,1)); % il y a des 1 l� ou j'ai des NaN
[row_nan,~] = find(TF_nan>0);
min_NaN = min(row_nan);
if size(min_NaN,1) == 0
    min_NaN = size(MatrixCentroid_t_brut,1);
end

MatrixCentroid_t_brut = MatrixCentroid_t_brut(1:min_NaN-1,:,:);

%% Create an image with the centroid on it

image_identifier_centroid = imread([pathMovie filesep 't' num2str(1,'%04d') filesep nameMovie '_t' num2str(1,'%04d') '_z' num2str(zMilieu_t(1,1)+1,'%04d') '.tif']);
image_identifier_centroid = imadjust(image_identifier_centroid);

% Initlialize mask des centroid
BW_centroid = zeros(size(image_identifier_centroid,1),size(image_identifier_centroid,2));

% Centroid = 1
for iteration = 1:size(MatrixCentroid_t_brut,1)
    
    if MatrixCentroid_t_brut(iteration,1,1)<1
        MatrixCentroid_t_brut(iteration,1,1) = 1;
    end
    
    if MatrixCentroid_t_brut(iteration,2,1)<1
        MatrixCentroid_t_brut(iteration,2,1) = 1;
    end
    
    %iteration
    BW_centroid(round(MatrixCentroid_t_brut(iteration,2,1)),round(MatrixCentroid_t_brut(iteration,1,1))) = 1;
end

% Dilate centroid
se = strel('disk',5);
BW_centroid = imdilate(BW_centroid,se);
BW_centroid = imcomplement(BW_centroid);
BW_centroid = uint16(BW_centroid);

% Mets les centroid sur l'image
image_identifier_centroid = image_identifier_centroid.*BW_centroid;
image_identifier_centroid = imadjust(image_identifier_centroid);

%% Click manually on the centroid of the cells and press enter
% if a file is already created, the step is automatically skipped
    
    exist_yes_no = exist([pathMovie filesep 'Data' filesep 'centroid_identite_cell' '.mat']);
    
    if exist_yes_no == 0
        figure
        imshow(image_identifier_centroid)
        [xChoix,yChoix] = ginput;
        xChoixScale = xChoix*scale1D; % mettre � la bonne �chelle
        yChoixScale = yChoix*scale1D;
        close
        
        % Utiliser les centroid pour sauver les images avec une identification
        position = [xChoix yChoix];
        value = 1:1:size(xChoix,1);%nombreCelluleEtudie;
        image_identification = insertText(image_identifier_centroid,position,value,'AnchorPoint','Center','FontSize',18);
        folder = [pathMovie filesep 'Data'];
        imwrite(image_identification,fullfile(folder,'identification_cellule.png'));
        
        % Save xChoixScale et yChoixScale
        save([pathMovie filesep 'Data' filesep 'centroid_identite_cell' '.mat'],'xChoixScale','yChoixScale','-v7.3');
        centroid_identite_cell = load([pathMovie filesep 'Data' filesep 'centroid_identite_cell' '.mat']);
        
    else % si le fichier existe, on le charge
        centroid_identite_cell = load([pathMovie filesep 'Data' filesep 'centroid_identite_cell' '.mat']);
    end

xChoixScale = centroid_identite_cell.xChoixScale;
yChoixScale = centroid_identite_cell.yChoixScale;
MatrixCentroid_click = [xChoixScale yChoixScale];

%% Scale of the images is adjusted to the real size

MatrixCentroid_t_brut = MatrixCentroid_t_brut*scale1D;
MatrixCentroid_t_brut_t1 = MatrixCentroid_t_brut(:,:,1);

%% Compare the click point and the object centroid in the image

matrice_R_click = NaN(size(MatrixCentroid_click,1),size(MatrixCentroid_click,2));

for iteration_cellule = 1:1:size(MatrixCentroid_t_brut_t1,1)
    
    xA2 = MatrixCentroid_click(:,1);
    yA2 = MatrixCentroid_click(:,2);
    
    xA1 = MatrixCentroid_t_brut_t1(iteration_cellule,1);
    yA1 = MatrixCentroid_t_brut_t1(iteration_cellule,2);
    
    r_carre_centroid = power((xA2-xA1),2)+power((yA2-yA1),2);
    r_centroid = power(r_carre_centroid,0.5);
    
    matrice_R_click(1:size(r_centroid,1),iteration_cellule) = r_centroid;
    
end

% Find the closest one
min_matrice_R_click = min(matrice_R_click);

% Select the n first with the smallest mistake
min_matrice_R_click = sort(min_matrice_R_click);
size_min = min(size(min_matrice_R_click,2),size(MatrixCentroid_click,1));
min_matrice_R_centroid_select = min_matrice_R_click(1:size_min);

%% Create final list of centroids

MatrixCentroid_final_t = NaN(size(MatrixCentroid_click,1),size(MatrixCentroid_click,2),size(MatrixCentroid_t_brut,3));

% Je retrouve ces valeurs dans la matrice
tf_t1 = ismember(matrice_R_click,min_matrice_R_centroid_select); % gives what is connected to what
[row_t1,col_t1] = find(tf_t1>0);

for iteration_temps = tmin:tmax
    for iteration_cellule = 1:1:size(col_t1,1)
        
        % x car seconde dimension = 1
        MatrixCentroid_final_t(row_t1(iteration_cellule),1,iteration_temps)= MatrixCentroid_t_brut(col_t1(iteration_cellule),1,iteration_temps);
        
        % y car seconde dimension = 2
        MatrixCentroid_final_t(row_t1(iteration_cellule),2,iteration_temps)= MatrixCentroid_t_brut(col_t1(iteration_cellule),2,iteration_temps);
    end
end

save([pathMovie filesep 'Data' filesep 'centroid_identite_cell_TEMPS' '.mat'],'MatrixCentroid_final_t','-v7.3');