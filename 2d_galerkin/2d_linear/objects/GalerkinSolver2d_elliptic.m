classdef (Abstract) GalerkinSolver2d_elliptic < GalerkinSolver2d

	properties
	end

	methods
		function self = GalerkinSolver2d_elliptic(dom,auxfun)
			
			%if nargin == 2
				% call superclass constructor
				self@GalerkinSolver2d(dom,auxfun);

				% calculate solution
				% ...
			%end

		end

		function self = solve(self)

			% assemble problem
			self = self.assembleTensors;
			self = self.assembleVectors;
			self = self.assembleBCs;
			[S,b] = self.finalAssembly;

			% load variables
			FreeNodes = self.domain.boundary.freeNodes;

			% solve and store solution
			v = sparse(self.domain.mesh.nNodes,1);
			v(FreeNodes) = S(FreeNodes,FreeNodes) \ b(FreeNodes);
			self.solution = v + self.vectors.U_D;

			% cleanup
			self = self.cleanup;

			%{
			% MANUAL TEST: STRATEGY 8

			% case 1: added by final assembly
			[S_1,b_1] = self.finalAssembly_test;

			% case 2: added here
			[S_2,b_2] = self.finalAssembly;
			%S_2 = S_2 + self.tensors.P;

			% solve and store solution
			v = sparse(self.domain.mesh.nNodes,1);
			v(FreeNodes) = S_2(FreeNodes,FreeNodes) \ b_2(FreeNodes);
			self.solution = v + self.vectors.U_D;

			% cleanup
			self = self.cleanup;
			%}
			

			%{
			% MANUAL TEST: STRATEGY 7
			% solve and store solution
			[S,b] = self.finalAssembly;
			S1 = S + self.tensors.P;
			[S2,b] = self.finalAssembly_test;


			FreeNodes = [7 8 10 11 12 13];
			v = sparse(self.domain.mesh.nNodes,1);
			v(FreeNodes) = S(FreeNodes,FreeNodes) \ b(FreeNodes);
			full(v(FreeNodes))
			%v(FreeNodes) = S2(FreeNodes,FreeNodes) \ b(FreeNodes);
			self.solution = v + self.vectors.U_D;

			S_test = S(FreeNodes,FreeNodes);
			full(S_test)
			S_test(1,1) = 2;
			S_test(1,5) = 0;
			S_test(1,6) = 0;
			full(S_test)
			v_test = S_test \ b(FreeNodes);
			full(v_test)

			full(S)

			% cleanup
			self = self.cleanup;
			%}

			%{
			% MANUAL TEST: STRATEGY 6
			S_correction_rows = sparse(zeros(size(S)));
			S_correction_cols = sparse(zeros(size(S)));
			S_correction_diag = sparse(zeros(size(S)));
			S_correction_rows(P_free,:) = S(P_rep,:);
			S_correction_cols(:,P_free) = S(:,P_rep);
			S_correction_diag(P_free,:) = S_correction_cols(P_rep,:);
			S_correction = S_correction_rows + S_correction_cols + S_correction_diag;
			%S = S + self.tensors.P;
			full(S)

			b_correction = zeros(size(b));
			b_correction(P_free) = b(P_rep);
			b = b + b_correction;

			% solve and store solution
			v = sparse(self.domain.mesh.nNodes,1);
			v(FreeNodes) = S(FreeNodes,FreeNodes) \ b(FreeNodes);
			self.solution = v + self.vectors.U_D;

			% cleanup
			self = self.cleanup;
			%}


			%{
			% MANUAL TEST: STRATEGY 5
			% In this test, the alterations to the stiffness matrix and source
			% vector are computed as separate tensor/vector and then added on to
			% the original.
			%FreeNodes = [8 10 11 12 13];
			%P_free = [1 7 2];
			%P_rep  = [5 9 4]

			%S_correction_rows = sparse(zeros(size(S)));
			%S_correction_cols = sparse(zeros(size(S)));
			%S_correction_rows(P_free,:) = S(P_rep,:);
			%S_correction_cols(:,P_free) = S(:,P_rep);
			%S_correction = S_correction_rows + S_correction_cols;
			%S = S + S_correction;

			%b_correction = zeros(size(b));
			%b_correction(P_free) = b(P_rep);
			%b = b + b_correction;

			%v = sparse(self.domain.mesh.nNodes,1);
			%v(FreeNodes) = S(FreeNodes,FreeNodes) \ b(FreeNodes);
			%self.solution = v + self.vectors.U_D;

			%self.solution(P_rep) = self.solution(P_free);
			%}
			
			%{
			% MANUAL TEST: STRATEGY 4
			% In this test, the periodic nodes are (manually) collected in
			% advance and the periodic matrix is made from alterations to the
			% original stiffness matrix.
			FreeNodes = [1 2 7 8 10 11 12 13];
			P_free = [1 7 2];
			P_rep  = [5 9 4]

			S_per = S;
			S_per(P_free,:) = S_per(P_free,:) + S_per(P_rep,:);
			S_per(:,P_free) = S_per(:,P_free) + S_per(:,P_rep);
			b(P_free) = b(P_free) + b(P_rep);

			v = sparse(self.domain.mesh.nNodes,1);
			v(FreeNodes) = S_per(FreeNodes,FreeNodes) \ b(FreeNodes);
			self.solution = v + self.vectors.U_D;

			self.solution(P_rep) = self.solution(P_free);
			%}

			%{
			% MANUAL TEST: STRATEGY 3
			% In this test the stiffness matrix for the problem is automatically generated, then the same alterations are carried out to make the stiffness matrix for the periodic problem.
			FreeNodes = [1 2 7 8 10 11 12 13];
			S_per = S;
			S_per(1,:) = S_per(1,:) + S_per(5,:);
			S_per(2,:) = S_per(2,:) + S_per(4,:);
			S_per(7,:) = S_per(7,:) + S_per(9,:);

			S_per(:,1) = S_per(:,1) + S_per(:,5);
			S_per(:,2) = S_per(:,2) + S_per(:,4);
			S_per(:,7) = S_per(:,7) + S_per(:,9);

			b(7) = b(7) + b(9);
			b(1) = b(1) + b(5);
			b(2) = b(2) + b(4);

			v = sparse(self.domain.mesh.nNodes,1);
			v(FreeNodes) = S_per(FreeNodes,FreeNodes) \ b(FreeNodes);
			self.solution = v + self.vectors.U_D;

			self.solution(5) = self.solution(1);
			self.solution(4) = self.solution(2);
			self.solution(9) = self.solution(7);
			%}

			%{
			% MANUAL TEST: STRATEGY 2
			% in this test, the stiffness matrix for the regular problem is
			% manually created, and then modified to create the stiffness matrix
			% for the periodic problem.
			FreeNodes = [1 2 7 8 10 11 12 13];
			S_reg = ...
			[[ 3  0  0  0  0 -1 -1  0  0  0 -1  0  0];
			 [ 0  3 -1  0  0  0 -1  0  0 -1  0  0  0];
			 [ 0 -1  5 -1  0  0  0 -1  0 -1  0 -1  0];
			 [ 0  0 -1  3  0  0  0  0 -1  0  0 -1  0];
			 [ 0  0  0  0  3 -1  0  0 -1  0  0  0 -1];
			 [-1  0  0  0 -1  5  0 -1  0  0 -1  0 -1];
			 [-1 -1  0  0  0  0  4  0  0 -1 -1  0  0];
			 [ 0  0 -1  0  0 -1  0  6  0 -1 -1 -1 -1];
			 [ 0  0  0 -1 -1  0  0  0  4  0  0 -1 -1];
			 [ 0 -1 -1  0  0  0 -1 -1  0  5 -1  0  0];
			 [-1  0  0  0  0 -1 -1 -1  0 -1  5  0  0];
			 [ 0  0 -1 -1  0  0  0 -1 -1  0  0  5 -1];
			 [ 0  0  0  0 -1 -1  0 -1 -1  0  0 -1  5]]; 

			S_per = S_reg;
			S_per(1,:) = S_per(1,:) + S_per(5,:);
			S_per(2,:) = S_per(2,:) + S_per(4,:);
			S_per(7,:) = S_per(7,:) + S_per(9,:);

			S_per(:,1) = S_per(:,1) + S_per(:,5);
			S_per(:,2) = S_per(:,2) + S_per(:,4);
			S_per(:,7) = S_per(:,7) + S_per(:,9);

			b(7) = b(7) + b(9)
			b(1) = b(1) + b(5)
			b(2) = b(2) + b(4)

			v = sparse(self.domain.mesh.nNodes,1);
			v(FreeNodes) = S_per(FreeNodes,FreeNodes) \ b(FreeNodes);
			self.solution = v + self.vectors.U_D;

			self.solution(5) = self.solution(1);
			self.solution(4) = self.solution(2);
			self.solution(9) = self.solution(7);
			%}

			%{
			% MANUAL TEST: STRATEGY 1
			% in this test, the stiffness matrix for periodic boundary
			% conditions is manually created. 
			FreeNodes = [1 2 7 8 10 11 12 13];
			FreeNodes = [8 10 11 12 13];
			S_reg = ...
			[[ 3  0  0  0  0 -1 -1  0  0  0 -1  0  0];
			 [ 0  3 -1  0  0  0 -1  0  0 -1  0  0  0];
			 [ 0 -1  5 -1  0  0  0 -1  0 -1  0 -1  0];
			 [ 0  0 -1  3  0  0  0  0 -1  0  0 -1  0];
			 [ 0  0  0  0  3 -1  0  0 -1  0  0  0 -1];
			 [-1  0  0  0 -1  5  0 -1  0  0 -1  0 -1];
			 [-1 -1  0  0  0  0  4  0  0 -1 -1  0  0];
			 [ 0  0 -1  0  0 -1  0  6  0 -1 -1 -1 -1];
			 [ 0  0  0 -1 -1  0  0  0  4  0  0 -1 -1];
			 [ 0 -1 -1  0  0  0 -1 -1  0  5 -1  0  0];
			 [-1  0  0  0  0 -1 -1 -1  0 -1  5  0  0];
			 [ 0  0 -1 -1  0  0  0 -1 -1  0  0  5 -1];
			 [ 0  0  0  0 -1 -1  0 -1 -1  0  0 -1  5]]; 

			S_per = ...
			[[ 6  0  0  0  3 -2 -1  0 -1  0 -1  0 -1];
			 [ 0  6 -2  3  0  0 -1  0 -1 -1  0 -1  0];
			 [ 0 -1  5 -1  0  0  0 -1  0 -1  0 -1  0];
			 [ 0  3 -2  3  0  0 -1  0 -1 -1  0 -1  0];
			 [ 0  0  0  0  6 -2 -1  0 -1  0 -1  0 -1];
			 [-1  0  0  0 -1  5  0 -1  0  0 -1  0 -1];
			 [-1 -1  0 -1 -1  0  8  0  8 -1 -1 -1 -1];
			 [ 0  0 -1  0  0 -1  0  6  0 -1 -1 -1 -1];
			 [-1 -1  0 -1 -1  0  4  0  4 -1 -1 -1 -1];
			 [ 0 -1 -1  0  0  0 -1 -1  0  5 -1  0  0];
			 [-1  0  0  0  0 -1 -1 -1  0 -1  5  0  0];
			 [ 0  0 -1 -1  0  0  0 -1 -1  0  0  5 -1];
			 [ 0  0  0  0 -1 -1  0 -1 -1  0  0 -1  5]]; 

			b(7) = b(7) + b(9)
			b(1) = b(1) + b(5)
			b(2) = b(2) + b(4)

			v = sparse(self.domain.mesh.nNodes,1);
			v(FreeNodes) = S_reg(FreeNodes,FreeNodes) \ b(FreeNodes);
			self.solution = v + self.vectors.U_D;

			self.solution(5) = self.solution(1);
			self.solution(4) = self.solution(2);
			self.solution(9) = self.solution(7);
			%}


		end

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

	end
end
