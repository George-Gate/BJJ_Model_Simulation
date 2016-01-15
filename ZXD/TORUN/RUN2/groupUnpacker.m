% The batch will call this function to run simulation for each parameter
% set
function groupUnpacker( pGroup, counter, outList )

    for i=1:counter
        results=adiabaticEvolution_Fock_reverse_RUN2(pGroup{i},outList);
        save(pGroup{i}.saveFile,'results');
    end

end

