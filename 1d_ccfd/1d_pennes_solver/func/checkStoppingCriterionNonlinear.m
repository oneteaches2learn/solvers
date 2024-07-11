function [U,Uiter,iterCase,iter] = checkStoppingCriterionNonlinear(...
					U,Uiter,dx,iterParm,i,iter)

	% unpack iteration parameters
	[atol,maxIter,verb] = constructIteration(iterParm);

	% compute residual
    residual = sqrt(sum(dx' .* (Uiter(2,:) - Uiter(1,:)).^2));

    % If convergence was acheived
    if residual < atol
		iterCase = 0;
		U(i,:)   = Uiter(2,:);
		if verb == 1
        	fprintf(' Timestep %g converged in %g iterations.\n',i-1,iter);
   		end
        
    % Else if too many iterations have occurred
    elseif iter >= maxIter 
		iterCase = 1;
        fprintf('   Timestep %g FAILED to converge. \n\n',i-1); 
        
    % Else prepare for next iteration    
    else
		iterCase = 2;
        iter = iter+1; 
        Uiter(1,:) = Uiter(2,:); 
    end

end
