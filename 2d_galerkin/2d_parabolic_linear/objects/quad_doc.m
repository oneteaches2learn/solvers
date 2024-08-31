		% THREE POINT QUADRATURE
		function int = threePointQuadrature(self,f)

			% if f is a function
			if isa(f,'function_handle') || isa(f,'symfun')

				% convert symfun to function_handle
				if isa(f,'symfun'), f = matlabFunction(f); end

				% get edge midpoints
				edges = self.getMeshEdges;

				% get edge midpoints
				midpts = self.getMeshEdgeMidpoints(edges);

				% compute F on edge midpoints 
				F = f(midpts(:,1),midpts(:,2));

			% else if f is a vector of midpoint values
			elseif isvector(f)

				F = f; 

				% convert to column, if necessary
				if isrow(F), F = F'; end 

			end

			% compute quadrature
			elemAvg = sum([F(self.elemEdges(:,1)), ...
					F(self.elemEdges(:,2)),F(self.elemEdges(:,3))],2) / 3;
			int = sum(self.elemAreas' .* elemAvg,"omitnan");

		end

		function IP = L2_IP_threePointQuadrature(self,arg1,arg2)

			% process inputs into vectors
			vals = Domain2d@function2midpoints(arg1,arg2);
			Arg1 = vals(:,1);
			Arg2 = vals(:,2);

			% Compute the inner product
			IP = self.threePointQuadrature(Arg1 .* Arg2);

		end	

		function norm = L2norm_threePointQuadrature(self,arg)

			% process input into vector
			Arg = Domain2d@function2midpoints(arg);

			% compute the norm
			norm_squared = self.L2_IP_threePointQuadrature(Arg,Arg);
			norm = sqrt(norm_squared);

		end

		function err = L2err_threePointQuadrature(self,arg1,arg2)

			% process inputs into vectors
			vals = Domain2d@function2midpoints(arg1,arg2);
			Arg1 = vals(:,1);
			Arg2 = vals(:,2);

			% compute error
			err = self.L2norm_threePointQuadrature(Arg1 - Arg2);

		end

		function vals = function2midpoints(self,varargin)

			temp = varargin;

			% convert symfun to function_handle
			for i = 1:length(temp)
				if isa(temp{i},'symfun')
					temp{i} = matlabFunction(temp{i});
				end
			end

			% compute function_handles on midpoints
			for i = 1:length(temp)
				if isa(temp{i},'function_handle')
					func = temp{i};
					temp{i} = func(self.mesh.midpoints);
				end
			end

			% check that values are column vectors
			for i = 1:length(temp)
				if isrow(temp{i})
					temp{i} = temp{i}';
				end
			end

			% store output as columns of matrix
			for i = 1:length(temp)
				vals(:,i) = temp{i};
			end

		end










		function IP = L2_IP_threePointQuadrature_nodal(self,arg1,arg2)

			[Arg1,Arg2] = self.nodal2midpoint(arg1,arg2);
			IP = self.L2_IP_threePointQuadrature(Arg1,Arg2);

		end

		function norm = L2norm_threePointQuadrature_nodal(self,f)

			% compute the norm
			norm_squared = self.L2_IP_threePointQuadrature_nodal(f,f);
			norm = sqrt(norm_squared);

		end


		function err = L2err_threePointQuadrature_nodal(self,arg1,arg2)

			[arg1,arg2,midpts] = self.nodal2midpoint(arg1,arg2);
			arg.arg1 = arg1;
			arg.arg2 = arg2;

			for i = 1:2
				if isa(arg(j),'symfun'), arg(j) = matlabFunction(arg(j)); end
			end

			if isa(arg(1),'function_handle') || isa(arg(2),'function_handle')
				edges = self.getMeshEdges;
				midpts = self.getMeshEdgeMidpoints(edges);
			end

			for i = 1:2
				if isa(arg(j),'function_handle')
					func = arg(j);
					arg(j) = func(midpts(:,1),midpts(:,2));
				end
			end

			err = self.L2norm_threePointQuadrature(arg.arg1 - arg.arg2);

		end


