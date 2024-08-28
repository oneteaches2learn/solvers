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

			elemIDs = setdiff([1:1:self.nElems],self.elements_Omega);

		end

		function nodeIDs = nodes_Omega(self)

			mesh = self.Mesh;
			nodeIDs_lower = findNodes(mesh,"region",Face=1);
			nodeIDs_upper = [];
			nodeIDs = unique([nodeIDs_lower,nodeIDs_upper]);

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
			mesh = self.Mesh;

			% plot all nodes
			h = pdemesh(mesh);

			% capture boundary data to redraw boundaries later
			boundaries = h(2);
			xBdry = boundaries.XData;
			yBdry = boundaries.YData;
			
			% having captured data, clear plot
			clf;

			% plot effective elements
			hold on
			h = pdemesh(mesh.Nodes,mesh.Elements(:,self.effectiveElems),...
						ElementLabels=x.ElementLabels);

			% Label nodes, optional
			if x.NodeLabels == "on"
				nNodes = size(self.Mesh.Nodes(self.effectiveNodes),2);	
				xData = self.Mesh.Nodes(1,self.effectiveNodes);
				yData = self.Mesh.Nodes(2,self.effectiveNodes);
				for i = 1:nNodes
					text(xData(i),yData(i),strcat('n',num2str(i)), ...
						'FontSize',x.NodeFontSize,'FontWeight','bold');
				end
			end

			% color inclusion elements
			incElem = self.elements_Qeps;
			k = pdemesh(mesh.Nodes,mesh.Elements(:,incElem));
			try, k.Color = [0.7 0.7 0.7]; end

			% replot boundary lines
			plot(xBdry,yBdry,"Color","red");
			hold off

		end

		% UTILITY FUNCTIONS
		function dom_on = get_ON(self)

			dom_on = self;
			dom_on.dl = self.dl(:,1:4);
			dom_on.edges = self.edges(1:4);
			dom_on.nEdges = length(dom_on.edges);

		end

		function dom_off = get_OFF(self)

			dom_off = self;
			dom_off.effectiveNodes = dom_off.nodes_Omega;
			dom_off.effectiveElems = dom_off.elements_Omega;
			dom_off.unusedNodes = dom_off.get_unusedNodes;
			dom_off.unusedElems = dom_off.get_unusedElems;
			dom_off.nNodes = length(dom_off.effectiveNodes);
			dom_off.nElems = length(dom_off.effectiveElems);
			dom_off.freeNodes = intersect(dom_off.freeNodes,dom_off.effectiveNodes);

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
