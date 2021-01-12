function [nameMovie , info_function] = identifier_film(k_movie)

% List of movie to analyse
% Possibility to add parameter
% nameMovie is the name of the movie (images has to be name accordingly)
% info_function is the file containing fundamental informations

if k_movie == 1
    nameMovie = 'movie1';
    info_function = str2func(['AIAinfo_movie1']);
    
elseif k_movie == 2
    nameMovie = 'movie2';
    info_function = str2func(['AIAinfo_movie2']);
end