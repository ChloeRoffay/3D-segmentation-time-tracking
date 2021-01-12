
function [MatrixCentroid_final_t_z] = Stage2_2_follow_cell_z(pathMovie,nameMovie,tmin,tmax,scale1D,nombreCelluleEtudie,zMilieu,k_animal,t,zmin,zmax,zMilieu_t);

%% Initialization

nombreCelluleMaxAttendues = 1000;
MatrixCentroid_brut = NaN(nombreCelluleMaxAttendues,2,zmax);
MatrixCentroid_z_brut = NaN(nombreCelluleMaxAttendues,2,zmax);

%% Extract the centroid of all the objects for every z

for z=zmin:zmax
    
    BW = imread([pathMovie filesep 't' num2str(t,'%04d') filesep 'Output_results' filesep 'Unionseg_' nameMovie '_t' num2str(t,'%04d') '_z' num2str(z,'%04d') '.png']);
    
    image_CC = bwconncomp(BW,4);
    Stat_All_Regions = regionprops(image_CC,'Centroid');
    
    Stat_All_Regions_cell = (struct2cell(Stat_All_Regions))';
    %matCentroid = cell2mat(Stat_All_Regions_cell(:,1))*scale1D; % VRAIE ECHELLE
    matCentroid = cell2mat(Stat_All_Regions_cell(:,1)); % % ECHELLE MATLAB
    
    MatrixCentroid_brut(1:size(matCentroid,1),1,z) = matCentroid(:,1);
    MatrixCentroid_brut(1:size(matCentroid,1),2,z) = matCentroid(:,2);
    
end

%% Extract distance of all the centroids compared to zMilieu + 2 

matrice_R_brut = NaN(nombreCelluleMaxAttendues,nombreCelluleMaxAttendues,tmax);

for iteration_z = zmin:zmax
    for iteration_cellule = 1:1:size(MatrixCentroid_brut,1)
        
        xA2 = MatrixCentroid_brut(:,1,iteration_z);
        yA2 = MatrixCentroid_brut(:,2,iteration_z);
        
        z_comparaison = zMilieu_t(t,1)+2;
        xA1 = MatrixCentroid_brut(iteration_cellule,1,z_comparaison);
        yA1 = MatrixCentroid_brut(iteration_cellule,2,z_comparaison);
        
        r_carre = power((xA2-xA1),2)+power((yA2-yA1),2);
        r = power(r_carre,0.5);
        
        matrice_R_brut(1:size(r,1),iteration_cellule,iteration_z) = r;
        
    end
end

%% Determine the closest centroid for each object of the image

% Select the minimum
min_matrice_R_brut = min(matrice_R_brut);

% Delete NaN
min_matrice_R_brut_tf = ~isnan(min_matrice_R_brut);
min_matrice_R_brut_tri = min_matrice_R_brut(min_matrice_R_brut_tf);

% Normal distribution of all the distances => mu is the max distance
% allowed to be considered as a candidate
pd_cellule = fitdist(min_matrice_R_brut_tri,'Normal');
r_max_autorise = 0.7*pd_cellule.sigma;
tf_r_autorise = min_matrice_R_brut_tri<r_max_autorise;

% Line number gives which centroid are connected (column number)
min_matrice_R_brut_tri = min_matrice_R_brut_tri(tf_r_autorise); 
tf_brut = ismember(matrice_R_brut,min_matrice_R_brut_tri);

%% Matrix is ordered accordingly

for iteration_z = zmin:zmax
    tf_temps = tf_brut(:,:,iteration_z);
    [row,col,~] = find(tf_temps>0);
    
    for iteration_cellule = 1:1:size(col,1)
        % x car seconde dimension = 1
        MatrixCentroid_z_brut(col(iteration_cellule),1,iteration_z)= MatrixCentroid_brut(row(iteration_cellule),1,iteration_z);
        
        % y car seconde dimension = 2
        MatrixCentroid_z_brut(col(iteration_cellule),2,iteration_z)= MatrixCentroid_brut(row(iteration_cellule),2,iteration_z);
    end
    
end

%% OSEF = delete the NaN
MatrixCentroid_z_brut_reshape = reshape(MatrixCentroid_z_brut,size(MatrixCentroid_z_brut,1),size(MatrixCentroid_z_brut,2)*size(MatrixCentroid_z_brut,3));

TF_nan = isnan(MatrixCentroid_z_brut_reshape); % il y a des 1 l� ou j'ai des NaN
[row_nan,~] = find(TF_nan<1); %je trouve les valeurs qui ne sont pas des NaN
min_NaN = max(row_nan);

if min_NaN<nombreCelluleEtudie
    min_NaN = nombreCelluleEtudie;
end

MatrixCentroid_z_brut = MatrixCentroid_z_brut(1:min_NaN,:,:);

%% Scale of the images is adjusted to the real size

MatrixCentroid_z_brut = MatrixCentroid_z_brut*scale1D;
MatrixCentroid_z_brut_zMilieu = MatrixCentroid_z_brut(:,:,zMilieu_t(t,1)+2);

%% Open the selected the centroid_t from before

fichier_temps = load([pathMovie filesep 'Data' filesep 'centroid_identite_cell_TEMPS' '.mat']);
liste_centroid_t_fixe = fichier_temps.MatrixCentroid_final_t(:,:,t);
MatrixCentroid_load = liste_centroid_t_fixe;

%% Compare the click point and the centroid in z

matrice_R_load = NaN(size(MatrixCentroid_load,1),size(MatrixCentroid_load,2));

for iteration_cellule = 1:1:size(MatrixCentroid_z_brut_zMilieu,1)
    
    xA2 = MatrixCentroid_load(:,1);
    yA2 = MatrixCentroid_load(:,2);
    
    xA1 = MatrixCentroid_z_brut_zMilieu(iteration_cellule,1);
    yA1 = MatrixCentroid_z_brut_zMilieu(iteration_cellule,2);
    
    r_carre_centroid = power((xA2-xA1),2)+power((yA2-yA1),2);
    r_centroid = power(r_carre_centroid,0.5);
    
    matrice_R_load(1:size(r_centroid,1),iteration_cellule) = r_centroid;
    
end

% Find the closest one
min_matrice_R_load = min(matrice_R_load,[],1);

% Select the n first with the smallest mistake
min_matrice_R_load = sort(min_matrice_R_load);
size_min = min(size(min_matrice_R_load,2),size(MatrixCentroid_load,1));
min_matrice_R_centroid_load = min_matrice_R_load(1:size_min);

%% Create final list of centroids

MatrixCentroid_final_t_z = NaN(size(MatrixCentroid_load,1),size(MatrixCentroid_load,2),size(MatrixCentroid_z_brut,3));

% Je retrouve ces valeurs dans la matrice
tf_t1 = ismember(matrice_R_load,min_matrice_R_centroid_load); % j'identifie ces valeurs => numero de ligne me donne quel centroid je veux connecter avec qui (numero de colonne)

% J'impose la taille pour �viter les probl�mes de taille dans le futur
row_t1 = NaN(size(tf_t1,1),1);
col_t1 = NaN(size(tf_t1,1),1);
[row_t1_inter,col_t1_inter] = find(tf_t1>0);
row_t1(1:size(row_t1_inter,1)) = row_t1_inter;
col_t1(1:size(col_t1_inter,1)) = col_t1_inter;

for iteration_z = zmin:zmax
    for iteration_cellule = 1:1:size(MatrixCentroid_load,1)
        
        if size(row_t1,1)==0
            continue
        end
        
        if isnan(row_t1(iteration_cellule)) == 1
            continue
        end
        
        % x car seconde dimension = 1
        MatrixCentroid_final_t_z(row_t1(iteration_cellule),1,iteration_z)= MatrixCentroid_z_brut(col_t1(iteration_cellule),1,iteration_z);
        
        % y car seconde dimension = 2
        MatrixCentroid_final_t_z(row_t1(iteration_cellule),2,iteration_z)= MatrixCentroid_z_brut(col_t1(iteration_cellule),2,iteration_z);
    end
end

