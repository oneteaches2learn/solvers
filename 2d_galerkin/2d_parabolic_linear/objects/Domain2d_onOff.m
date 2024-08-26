classdef Domain2d_onOff < Domain2d_punctured
% Domain2d_onOff is a domain with inclusions that can be turned on or off
%	
% Author: Tyler Fara					Date: August 21, 2024
%-----------------------------------------------------------------------------%

	methods
		% CONSTRUCTOR
		function self = Domain2d_onOff(x,y,inc,eps)
		
			% call Domain2d superclass constructor
			self@Domain2d_punctured(x,y,inc,eps);

			% create decomposed geometry description matrix 
			dl_domain = Domain2d_punctured.dl_domain(x,y);
			if nargin == 2
				dl_Qeps = [];
			else
				dl_Qeps = Domain2d_punctured.dl_Qeps(x,y,inc,eps);
			end
			dl = [dl_domain dl_Qeps];

			% modify dl to include Q^epsilon
			dl = Domain2d_onOff.modify_dl(dl);

			% store geometry matrix
			self.dl = dl;

		end

		% GETTERS
		function elemIDs = elements_Omega(self)
			
			mesh = self.Mesh;
			elemIDs = findElements(mesh,"region",Face=1);

		end

		function elemIDs = elements_Qeps(self)

			elemIDs = setdiff([1:1:self.nElem],self.elements_Omega);

		end

		% PLOTTERS
		function h = plot(self,NameValueArgs)

			arguments
				self
				NameValueArgs.NodeLabels string = "off"
				NameValueArgs.NodeFontSize double = 18
				NameValueArgs.ElementLabels string = "off"
			end
			x = NameValueArgs;

			% Plot mesh
			hold on
			h = pdemesh(self.Mesh,...
					ElementLabels=x.ElementLabels);

			% Label nodes, optional
			if x.NodeLabels == "on"
				nNodes = size(self.Mesh.Nodes,2);	
				xData = self.Mesh.Nodes(1,:);
				yData = self.Mesh.Nodes(2,:);
				for i = 1:nNodes
					text(xData(i),yData(i),strcat('n',num2str(i)), ...
						'FontSize',x.NodeFontSize,'FontWeight','bold');
				end
			end

			% capture boundary data to redraw boundaries later
			boundaries = h(2);
			xBdry = boundaries.XData;
			yBdry = boundaries.YData;

			% color inclusion elements
			mesh = self.Mesh;
			incElem = self.elements_Qeps;
			k = pdemesh(mesh.Nodes,mesh.Elements(:,incElem));
			k.Color = [0.7 0.7 0.7];

			% replot boundary lines
			plot(xBdry,yBdry,"Color","red");
			hold off

		end

	end

	methods (Static)
		function dl = modify_dl(dl)

			% get number of inclusions
			nInc = size(dl,2) / 4 - 1;

			incLabels = [1:1:nInc] + 1;
			incLabels = repmat(incLabels,4,1);
			incLabels = reshape(incLabels,1,[]);

			dl(6,5:end) = incLabels;

		end
	end
end
