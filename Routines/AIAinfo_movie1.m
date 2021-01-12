function   [t_choc,zmin,zmax,zMilieu,tmin,tmax,pathMovie,scale1D,scale2D,nombreCelluleEtudie,confiance,step_time] = AIAinfo_movie1(animalChloe)

% Path of the images
pathMovie = ['Z:\Chloe\Use\Traitement\20191126-latA-50nM-1h-hyper\' animalChloe];

zmin = 1;
zmax = 40;

tmin = 1;
tmax = 10;

t_choc = 1;

% Number of cell in the stack
nombreCelluleEtudie = 10;

step_time = 60; % in seconds

% zMilieu = medium plan of the cell
% Plan where you see all the cells
zMilieu = 7;

% Pixel size in microns
binning = 4;
scale1D = 0.065*binning;
scale2D = scale1D^2;

confiance = 5;
% Confiance en ma capacit� � identifier le centroid d'une cellule
% Augmenter la confiance si je ne trouve pas toutes les cellules ou s'il manque des bouts

