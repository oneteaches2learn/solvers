classdef NewtonRaphsonSolver2d_poisson < NewtonRaphsonSolver2d_elliptic

	properties
	end

	methods
		function self = NewtonRaphsonSolver2d_poisson(dom,auxfun)
			
			% call superclass constructor
			self@NewtonRaphsonSolver2d_elliptic(dom,auxfun);

		end

		function self = assembleTensors(self,U)

			% Initialisation
			nodes = self.domain.mesh.nodes;
			elems = self.domain.mesh.elements;

			A = sparse(self.domain.mesh.nNodes,self.domain.mesh.nNodes);
			for j = 1:self.domain.mesh.nElems
				A(elems(j,:),elems(j,:)) = A(elems(j,:),elems(j,:)) ...
					+ self.localdj(nodes(elems(j,:),:),U(elems(j,:)));
			end

			self.tensors.DJ = A;

		end

		function self = assembleVectors(self,U)

			% Initialisation
			nodes = self.domain.mesh.nodes;
			elems = self.domain.mesh.elements;
			[temp,b_cubic] = self.domain.nodalQuadrature(U.^3);

			b = sparse(size(nodes,1),1);
			for j = 1:self.domain.mesh.nElems
				b(elems(j,:)) = b(elems(j,:)) ...
					+ self.localj(nodes(elems(j,:),:),U(elems(j,:)),b_cubic(j));
			end

			self.vectors.b_diffeq = b;
			self.vectors.b_load = self.computeVolumeForces(U);

		end

		function self = assembleBCs(self)

			%self = self.computePeriodicBCs;
			self.vectors.U_D   = self.computeDirichletBCs;
			self.vectors.b_neu = self.computeNeumannBCs;
			[self.tensors.E,self.vectors.b_rob] = self.computeRobinBCs;

		end

		function [DJ,J] = finalAssembly(self,U)

			DJ = self.tensors.DJ + self.tensors.E;
			J = self.vectors.b_diffeq + self.tensors.E * U - self.vectors.b_load + self.vectors.b_neu - self.vectors.b_rob;

		end

		function M = localdj(self,vertices,U_loc)

			% compute local contributions to DJ
			A_loc = self.localStiffnessMatrix(vertices);
			Mr_loc_old = -self.localMassMatrix + self.localDCubic_direct(U_loc);
			Mr_loc_nodal = self.localQuadrature_nodal(self.coefficients.dr_du(vertices(:,1),vertices(:,2),U_loc));
			Mr_loc_centroid = self.localQuadrature_centroid(vertices,U_loc,self.coefficients.dr_du);
			Mr_loc = Mr_loc_old;

			% assemble local contributions to DJ
			M = A_loc + Mr_loc;
			Area = det([ones(1,3);vertices']) / 2;
			M = M * Area;
			
		end
		
		function j_loc = localj(self,vertices,U_loc,U_cubed)

			% compute local contributions to J
			A_loc = self.localStiffnessMatrix(vertices);
			b_r_old = -self.localMassMatrix * U_loc + self.localCubic_direct(U_loc);
			b_r = self.localQuadrature(self.coefficients.r(vertices(:,1),vertices(:,2),U_loc));
			b_r = b_r_old;

			% assemble local contributions to J
			Area = det([ones(1,3); vertices']) / 2;
			j_loc = A_loc * U_loc + b_r;
			j_loc = Area * j_loc;

		end

		function b = localQuadrature(self,U)

			M_loc = [2,1,1; 1,2,1; 1,1,2] / 12;
			b = M_loc * U;
			
		end

		function A_loc = localStiffnessMatrix(self,vertices)

			% compute local stiffness matrix
			G = [ones(1,3);vertices'] \ [zeros(1,2);eye(2)];
			A_loc = G * G';

		end

		function M_r = localQuadrature_nodal(self,C_U)

			CU_nodal = sum(C_U) / 3;
			M_loc = [2,1,1; 1,2,1; 1,1,2] / 12;
			M_r = CU_nodal * M_loc;

		end

		function M_r = localQuadrature_centroid(self,X,U,c)

			X_centroid = sum(X) / 3;
			U_centroid = sum(U) / 3;
			CU_centroid = c(X_centroid(:,1),X_centroid(:,2),U_centroid);
			M_loc = [2,1,1; 1,2,1; 1,1,2] / 12;
			M_r = CU_centroid * M_loc;

		end

	end
end
