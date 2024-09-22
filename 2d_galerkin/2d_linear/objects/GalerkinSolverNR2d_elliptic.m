classdef GalerkinSolverNR2d_elliptic < GalerkinSolver2d

	properties
	end

	methods
		function self = GalerkinSolverNR2d_elliptic(dom,auxfun)
			
			% call superclass constructor
			self@GalerkinSolver2d(dom,auxfun);

			if nargin == 2

				% calculate solution
            	self = self.solve;

			end

		end

		function self = solve(self)
						
			% Initialisation
			coordinates = self.domain.mesh.nodes;
			elements3 = self.domain.mesh.elements;
			dirichlet = self.domain.boundary.D_nodes;
			neumann = [];
			FreeNodes = self.domain.boundary.freeNodes;

			% Initial value
			U = -ones(size(coordinates,1),1);
			U(unique(dirichlet)) = u_d(coordinates(unique(dirichlet),:));

			% Newton-Raphson iteration
			for i = 1:50
			
				% Assembly of DJ(U)
				A = sparse(size(coordinates,1),size(coordinates,1));
				for j = 1:size(elements3,1)
					A(elements3(j,:),elements3(j,:)) = A(elements3(j,:),elements3(j,:)) ...
					+ self.localdj(coordinates(elements3(j,:),:),U(elements3(j,:)));
				end
				
				% Assembly of J(U)
				b = sparse(size(coordinates,1),1);
				for j = 1:size(elements3,1)
					b(elements3(j,:)) = b(elements3(j,:)) ...
					+ self.localj(coordinates(elements3(j,:),:),U(elements3(j,:)));
				end
				
				% Volume Forces
				for j = 1:size(elements3,1)
					b(elements3(j,:)) = b(elements3(j,:)) + ...
					det([1 1 1; coordinates(elements3(j,:),:)']) * ...
					f(sum(coordinates(elements3(j,:),:))/3)/6;
				end
				
				% Neumann conditions
				for j = 1 : size(neumann,1)
					b(neumann(j,:))=b(neumann(j,:)) - norm(coordinates(neumann(j,1),:)- ...
					coordinates(neumann(j,2),:))*g(sum(coordinates(neumann(j,:),:))/2)/2;
				end
				
				% Dirichlet conditions
				W = zeros(size(coordinates,1),1);
				W(unique(dirichlet)) = 0;
				
				% Solving one Newton step
				W(FreeNodes) = A(FreeNodes,FreeNodes)\b(FreeNodes);
				U = U - W;
				if norm(W) < 10^(-10)
					break
				end
			end

			% graphic representation
			% show(elements3,[],coordinates,full(U));
			self.solution = U;

		end

		function out = testfunc(self,x)

			out = sin(2 * pi * x(1)) * sin(2 * pi * x(2));

		end

		function M = localdj(self,vertices,U)

			Eps = 1/100;
			G = [ones(1,3);vertices'] \ [zeros(1,2);eye(2)];
			Area = det([ones(1,3);vertices']) / 2;
			M = Area*(Eps*G*G'-[2,1,1;1,2,1;1,1,2]/12 + ...
				[12*U(1)^2+2*(U(2)^2+U(3)^2+U(2)*U(3))+6*U(1)*(U(2)+U(3)),...
				3*(U(1)^2+U(2)^2)+U(3)^2+4*U(1)*U(2)+2*U(3)*(U(1)+U(2)),...
				3*(U(1)^2+U(3)^2)+U(2)^2+4*U(1)*U(3)+2*U(2)*(U(1)+U(3));
			3*(U(1)^2+U(2)^2)+U(3)^2+4*U(1)*U(2)+2*U(3)*(U(1)+U(2)),...
				12*U(2)^2+2*(U(1)^2+U(3)^2+U(1)*U(3))+6*U(2)*(U(1)+U(3)),...
				3*(U(2)^2+U(3)^2)+U(1)^2+4*U(2)*U(3)+2*U(1)*(U(2)+U(3));
			3*(U(1)^2+U(3)^2)+U(2)^2+4*U(1)*U(3)+2*U(2)*(U(1)+U(3)),...
				3*(U(2)^2+U(3)^2)+U(1)^2+4*U(2)*U(3)+2*U(1)*(U(2)+U(3)),...
				12*U(3)^2+2*(U(1)^2+U(2)^2+U(1)*U(2))+6*U(3)*(U(1)+U(2))]/60);

		end

		function b = localj(self,vertices,U)

			Eps = 1/100;
			G = [ones(1,3);vertices'] \ [zeros(1,2);eye(2)];
			Area = det([ones(1,3);vertices']) / 2;
			b = Area*((Eps*G*G'-[2,1,1;1,2,1;1,1,2]/12)*U+ ...
				[4*U(1)^3+ U(2)^3+U(3)^3+3*U(1)^2*(U(2)+U(3))+2*U(1) ...
				*(U(2)^2+U(3)^2)+U(2)*U(3)*(U(2)+U(3))+2*U(1)*U(2)*U(3);
			4*U(2)^3+ U(1)^3+U(3)^3+3*U(2)^2*(U(1)+U(3))+2*U(2) ...
				*(U(1)^2+U(3)^2)+U(1)*U(3)*(U(1)+U(3))+2*U(1)*U(2)*U(3);
			4*U(3)^3+ U(2)^3+U(1)^3+3*U(3)^2*(U(2)+U(1))+2*U(3) ...
				*(U(2)^2+U(1)^2)+U(2)*U(1)*(U(2)+U(1))+2*U(1)*U(2)*U(3)]/60);

		end
		

		function VolumeForce = source(self,x);

			x
			VolumeForce = zeros(size(x,1),1);

		end
		

		function Stress = g(self,x)

			Stress = zeros(size(x,1),1);

		end
		
		%{
		function self = assembleTensors(self) 

			self.tensors.A = self.assembleStiffnessMatrix;
			self.tensors.M_r = self.assembleMassMatrix(self.coefficients.r);

		end

		function self = assembleVectors(self)

			self.vectors.b_vol = self.computeVolumeForces();

		end

		function [S,b] = finalAssembly(self)

			% NOTE: placeholder function. Actually handled by specific subclasses.
			...

		end
		%}

	end
end
