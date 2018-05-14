% SCRIPT: makeEdgeletLib
%    - creates the .mat file and variable mpLibrary by calling makeEdgelet repeatedly 
% mpLibrary is a sparse 4-D cell array indexed like:
%		mpLibrary{type, size, orientation, version}
%		where type = 1 refers to narrowband
%		      type = 2 refers to broadband
%		      size & ori are actual values 
%		eg: mpLibrary{1, 64, 360, 1} is a valid call
%           mpLibrary{1, 64, 360, 16} is a valid call only if phase or
%           orientation are scrambled
% change these to include more sizes or orientations in the library
% size should not be smaller than 16 due to resolution issues

s = [4 8 16 32];
o = [0:30:330];

nType = 9;
nSize = length(s); 
nOri  = length(o);
nVari = 25;

mpLibrary = cell(nType, nSize, nOri);

for i = 1:nOri   
    for j = 1:nSize %generate one variation if there is no randomness
        mpLibrary{1, j, i, 1} = zeros(s(j));
        mpLibrary{1, j, i, 1} = makeEdgelet_os('narrowband', s(j), o(i),0,0);
        
        mpLibrary{2, j, i, 1} = zeros(s(j));
        mpLibrary{2, j, i, 1} = makeEdgelet_os('broadband', s(j), o(i),0,0);
            
        mpLibrary{3, j, i, 1} = zeros(s(j));
        mpLibrary{3, j, i, 1} = makeEdgelet_os('gabor', s(j), o(i),0,0);
        for k = 1:nVari %generate nVari variations if there is randomness
            mpLibrary{4, j, i, k} = zeros(s(j));
            mpLibrary{4, j, i, k} = makeEdgelet_os('narrowband', s(j), o(i),1,0);
        
            mpLibrary{5, j, i, k} = zeros(s(j));
            mpLibrary{5, j, i, k} = makeEdgelet_os('broadband', s(j), o(i),1,0);
            
            mpLibrary{6, j, i, k} = zeros(s(j));
            mpLibrary{6, j, i, k} = makeEdgelet_os('gabor', s(j), o(i),1,0);
            
            mpLibrary{7, j, i, k} = zeros(s(j));
            mpLibrary{7, j, i, k} = makeEdgelet_os('narrowband', s(j), o(i),0,1);
        
            mpLibrary{8, j, i, k} = zeros(s(j));
            mpLibrary{8, j, i, k} = makeEdgelet_os('broadband', s(j), o(i),0,1);
            
            mpLibrary{9, j, i, k} = zeros(s(j));
            mpLibrary{9, j, i, k} = makeEdgelet_os('gabor', s(j), o(i),0,1);
        end
    end
end

save 'mpLibrary_os.mat' mpLibrary s o nVari;