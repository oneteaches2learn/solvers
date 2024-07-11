function [atol,maxIter,verb] = constructIteration(iterParm);

	if isa(iterParm,'char') && (iterParm == 'default');
		atol    = 1e-6;
		maxIter = 20;
		verb    = 0; 
	else
		atol    = iterParm{1}; 
		maxIter = iterParm{2}; 
		verb    = iterParm{3}; 
	end

end


