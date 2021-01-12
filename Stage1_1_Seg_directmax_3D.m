function Stage1_1_Seg_directmax_3D(t,nameMovie,pathMovie,zmin,zmax,tmin,zMilieu,nombreCelluleEtudie)
%% Parameter to enter

% Intensity of blurring
% Higher => more blurry => undersegmented
blur_intensity = 5;

% Manual choices required
% verif_segmentation = 1 to check the segmentation quality on one slice
% verif_segmentation = 0 to not do anything
verif_segmentation = 0;

%% Parameter needed

nombre_z = zmax-zmin+1;

% connectivity for the research of region of minima
connect_extmin= 26; % 8 in 2D, 26 in 3D

% connectivity for the algorithm of watershed
connect_watershed= 26; % 8 in 2D, 26 in 3D

% Create the file Output_results for each t
dossier_save_dirsegimage = [pathMovie filesep 't' num2str(t,'%04d') filesep 'Output_results' ];
status = mkdir(dossier_save_dirsegimage);

dossier_save_keyframe = [pathMovie filesep 't' num2str(t,'%04d') filesep 'Output_backup' ];
status = mkdir(dossier_save_keyframe);

%% FAIRE LE STACK

% Read size of the image
I_size= imread([pathMovie filesep 't' num2str(t,'%04d') filesep nameMovie '_t' num2str(t,'%04d') '_z' num2str(zmin,'%04d') '.tif']);

% Create matrix of the z-stack, fixed t
I_3D = zeros(size(I_size,1),size(I_size,2),nombre_z);
I_3D = uint16(I_3D);

% Blur the images
for z = zmin:1:zmax
    I = imread([pathMovie filesep 't' num2str(t,'%04d') filesep nameMovie '_t' num2str(t,'%04d') '_z' num2str(z,'%04d') '.tif']);
    I = imgaussfilt(I,blur_intensity);
    I_3D(:,:,z-zmin+1) = I;
end

% Adjust intensity on the full stack
% Smarter than adjusting slice by slice
I_3D = imadjust_stack(I_3D);

%% Find automatically segmentation threshold
% Depends on image quality, the probe, the microscope.......
% You have to play around at first

% Distribution of pixel intensities of the full stack
pd_cellule = fitdist(I_3D(:),'Normal');

% Cut at an empiric value
parameter_threshold = 0.2*pd_cellule.mu;

%% Segmentation step
% Everything is done on a full z-stack
% Longer processing but more consistent as well

% Search the region of extended minima
extmin= imextendedmin(I_3D,parameter_threshold,connect_extmin);
disp('extmin done')

% Creates an image that forces the previous regions to be minima
I2 = imimposemin(I_3D,extmin);
disp('imimposemin done')

% Simple watershed on the image with seeded minima
dirsegimage = ~watershed(I2,connect_watershed);
disp('watershed done')

%OPTION
dirsegimage = bwmorph3(~dirsegimage,'fill');

%% Check segmentation
answer = 1;
if verif_segmentation ==  1
        figure('units','normalized','outerposition',[0 0 1 1])
        
        % Plot the microscopy image
        subplot(1,3,1)
        just_look = imread([pathMovie filesep 't' num2str(t,'%04d') filesep nameMovie '_t' num2str(t,'%04d') '_z' num2str(zMilieu+2,'%04d') '.tif']);
        just_look = imadjust(just_look);
        imshow(just_look)
        
        % Plot the microscopy image + blurring
        subplot(1,3,2)
        just_look_blur = I_3D(:,:,zMilieu-zmin+1);
        imshow(just_look_blur)
        
        % Plot the skeleton
        subplot(1,3,3)
        imshow(dirsegimage(:,:,zMilieu-zmin+1))
        
        % Choose to continue or to stop
        choice = questdlg('Content de cette segmentation ?', ...
            'Dialog box', ...
            'Yes it is perfect','No, it is bad','No, it is bad');
        % Handle response
        switch choice
            case 'Yes it is perfect'
                disp([choice ' coming right up.'])
                answer = 1;
            case 'No, it is bad'
                disp([choice ' coming right up.'])
                answer = 0;
        end
end
%% Save the segmentation skeleton

% If segmentation is good, skeleton are saved
if answer == 0
    disp(['No segmentation... TRY AGAIN']);
else
    for z = zmin:1:zmax
        disp(['Creating Directskel file for frame #' num2str(t) ' and slice #' num2str(z) '...']);
        imwrite(dirsegimage(:,:,z-zmin+1),[pathMovie filesep 't' num2str(t,'%04d') filesep 'Output_results' filesep 'directskel_' nameMovie '_t' num2str(t,'%04d') '_z' num2str(z,'%04d') '.png'],'png');
    end
end

