%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Author: Chance and Charly
%
%  Function Discription : 
%  Given the number of exemplars desired (exp) the function will output
%  every possible texture an (exp) amount of time and add them to an array
%  (stimArray)to be used later. A sentinal value is also passed that checks
%  to see if the user only wants  first x amount of exemplar's.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Function: stimGen
%       Input:  
%               exp :   This is the number of exempars that the user want
%                       when genterating the textures
%               all :   This acts as a sentinal value: if all = 0 then the
%                       function will return all of the desired exemplars
%                       of every combinations of micorpatters, orientation,
%                       and size; if all > 0 the function will return the
%                       first 'all' textures. 
%       Output:
%               stimArray :     As the function loops through the for loop 
%                               it will systematically create every 
%                               possible texture combination, then stores 
%                               these textures in an array that is then 
%                               returned to the function caller.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [stimArray] = stimGen( exp , all )

load mpLibrary_os;
global mpLibrary s o nVari;
patterns = [];

x = 1;
for nPatt = 9 : 11
    for Oris = 0 : 11
        for size = 3 : 5
            for exemplars = 1 : exp
                fprintf('%d | [%d : %d : %d] \n',x, nPatt, Oris, size);
                myFirstTexture = TextureGenCH(64, 2^nPatt, 3, [Oris*30], 2^size);
                patterns{x} = myFirstTexture;
                
                x = x + 1;
                if (x > all & not(eq(all,0)))
                    stimArray = patterns;
                    return;
                end  
                
            end
        end    
    end
end
 stimArray = patterns;
 return;

end

