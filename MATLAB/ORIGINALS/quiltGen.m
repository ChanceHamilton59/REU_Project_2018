function [ quiltArray ] = quiltGen( stimArray, exemplars )

stim.modType   = 'quilt';
stim.modDepth  = 100; 
stim.carrA.con = 0.14; 
stim.carrB.con = 0.14;


Array = [];
index = 1; 
for cA = 1 : length(stimArray)
    stim.carrPattA = stimArray{cA};
       
    for cB = 1 : length(stimArray)
        if (cB ~= cA)
           stim.carrPattB = stimArray{cB};
           quilt = stimQuilter(stim, 0);
           Array{index} = quilt;
           %imagesc(Array{index}); scrollsubplot(10,5,index); colormap('gray');
           index = index + 1;
        end    

    end
end

for cA = 1 : length(stimArray)
    stim.carrPattA = stimArray{cA};
       
    for cB = 1 : length(stimArray)
        if (cB ~= cA)
           stim.carrPattB = stimArray{cB};
           quilt = stimQuilter(stim, 1);
           Array{index} = quilt;
           %imagesc(Array{index}); scrollsubplot(10,5,index); colormap('gray');
           index = index + 1;
        end    

    end
end
quiltArray = Array;

