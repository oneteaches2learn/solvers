classdef PhysiologicalBC
% PHYSIOLOGICALBC is a class for handling physiological boundary conditions
% in 2D Galerkin linear problems.
%
% This class assembles the boundary condition function from the provided
% physical data and then outputs the function and its derivative, which can be
% stored in the main problem class.
%
% MODEL FORMAT:
% The physiological boundary condition is of the form:
%
%    alpha(x,t,u,v) = conv(u) + rad(u) + evap(u,v)
%
% where:
%
%      conv(u) = c_1 * (u - u_inf)
%       rad(u) = c_2 * ((u + 273.15)^4 - (u_inf + 273.15)^4)
%    evap(u,v) = [c_3 * (u - u_*) + c_4 * (v - v_*)]_pos + c_5
%
% Here, c_1, c_2, c_3, c_4, c_5 are coefficients (which may be functions of
% space and time), u_inf is the ambient (i.e. "far-field") temperature, u_* and
% v_* are set point temperatures that the body tries to maintain. The notation
% [x]_pos is the positive part of x, i.e., defined as [x]_pos = max(x, 0). All
% temperatures are given in Celsius.
% 
%
% COEFFICIENTS:
% Per [1], coefficient c_1 is defined:
%
%     c_1 = B v^m
%
% where B and m are emperical constants (depending on body segment and body
% position, i.e. standing versus sitting) and v is the air velocity. In [1], for
% the human forearm, B = 12.5 and m = 0.54.
% 
% Per [1], coefficient c_2 is defined:
%
%     c_2 = view_factor * emissitivity * stefan_boltzmann
%
% where the view_factor = A_r / A_D with A_r being the effective radiative area
% of the body and A_D the DuBois area; this view_factor is given in [1] as 0.73
% for a standing human and 0.70 for a seated human. Emissivity is given in [1]
% as 0.95 for human skin. 
%
% Per [2], coefficients c_3 and c_4 are defined as:
%
%     c_3 = 74.4 * a_SW
%     c_4 = 31.9 * a_SW
% 
% where a_SW is a distribution factor for sweat rate and depends on the body
% segment. For the hand it is given in [2] as a_SW = 0.055. The coefficients
% 74.4 and 31.9 account for the relative importance of local versus global
% temperature errors compared to set points, i.e. in u - u_* and v - v_*,
% respectively. (Note that in [2], there is a Hill function dependence on the
% rate of change of skin temperature as well. That is, rapid heating of the skin
% can radically increase sweat rate. However, this term is set to zero for
% cooling skin, reflecting the physical impossibility of negative sweat rates.
% Since we assume a cooling environment, we neglect this term here.)
% 
% Per [2], coefficient c_5 is defined as follows:
%
%     c_5 = 5.6, if u_inf < 15, 
%         = 5.6 + 0.0213 * (u_inf - 15.0)^2, if 15 < u_inf < 30,
%         = 10.4, if u_inf > 30.
%
% For our application, we assume that u_inf is cold, so we set c_5 = 5.6.
%
%
% REFERENCES:
% [1] de Dear, Richard J.; Arens, Edward; Hui, Zhang; Oguro, Masayuki.
% "Convective and radiative heat transfer coefficients for individual human body
% segments", International Journal of Biometeorology, 1997. 
% [2] Castellani, Michael P.; Rioux, Timothy P.; Castellani, John W.; Potter,
% Adam W.; Xu, Xiaojiang. "A geometrically accurate 3 dimensional model of human
% thermoregulation for transient cold and hot environments." Computers in
% Biology and Medicine, 2021.

    properties
        view_factor = 0.73;             % standing human
        emissivity = 0.95;              % human skin
        stefan_boltzmann = 5.67e-8;     % W/m
        B = 12.5;                       % forearm, standing
        m = 0.54;                       % forearm, standing
        a_SW = 0.055;                   % hand
        uStar                           % set point skin temperature, Celsius
        vStar                           % set point core temperature, Celsius
        uInf                            % ambient temperature, Celsius
        windSpeed                       % wind speed, m/s
        coefficients                    % vector to hold coefficients
    end

    methods
        function self = PhysiologicalBC(uStar, vStar, uInf, windSpeed)

			arguments
                uStar double = 33.0
                vStar double = 37.0
                uInf double  = -40.0
                windSpeed double = 1.0
            end

            % store inputs
            self.uStar = uStar;
            self.vStar = vStar;
            self.uInf  = uInf;
            self.windSpeed = windSpeed;

        end
        
        function func = BC(self)

            % store individual components
            conv = self.convection();
            rad  = self.radiation();
            evap = self.evaporation();

            % assemble full BC
            func = @(x1,x2,t,u,v) (conv(x1,x2,t,u,v) + rad(x1,x2,t,u,v) + evap(x1,x2,t,u,v) );

        end

        function func = BC_du(self)

            % store indvidiual components
            conv = self.convection_du();
            rad  = self.radiation_du();
            evap = self.evaporation_du();

            % assemble full BC
            func = @(x1,x2,t,u,v) (conv(x1,x2,t,u,v) + rad(x1,x2,t,u,v) + evap(x1,x2,t,u,v) );

        end

        function func = BC_linearized(self)

            % store individual components
            conv = self.convection();
            rad  = self.radiation_linearized();
            evap = self.evaporation();

            % assemble full BC
            func = @(x1,x2,t,u,v) (conv(x1,x2,t,u,v) + rad(x1,x2,t,u,v) + evap(x1,x2,t,u,v) );

        end

        function func = BC_linearized_du(self)

            % store individual components
            conv = self.convection_du();
            rad  = self.radiation_linearized_du();
            evap = self.evaporation_du();

            % assemble full BC
            func = @(x1,x2,t,u,v) (conv(x1,x2,t,u,v) + rad(x1,x2,t,u,v) + evap(x1,x2,t,u,v) );

        end

        function func = convection(self)
            cofs = self.coefficients;
            func = @(x1,x2,t,u,v) cofs(1) .* (u - self.uInf);
        end

        function func = convection_du(self)
            cofs = self.coefficients;
            func = @(x1,x2,t,u,v) cofs(1);
        end

        function func = radiation(self)
            cofs = self.coefficients;
            func = @(x1,x2,t,u,v) cofs(2) .* ((u + 273.15).^4 - (self.uInf + 273.15).^4);
        end

        function func = radiation_du(self)
            cofs = self.coefficients;
            func = @(x1,x2,t,u,v) cofs(2) .* 4 .* (u + 273.15).^3;
        end

        function func = radiation_linearized(self)

            func = @(x1,x2,t,u,v) 4.27 .* (u - self.uInf);

        end

        function func = radiation_linearized_du(self)

            func = @(x1,x2,t,u,v) 4.27;

        end

        function func = evaporation(self)
            cofs = self.coefficients;
            func = @(x1,x2,t,u,v) max(cofs(3) .* (u - self.uStar) + cofs(4) .* (v - self.vStar),0) + cofs(5);
        end

        function func = evaporation_du(self)
            cofs = self.coefficients;
            func = @(x1,x2,t,u,v) (cofs(3) .* (u - self.uStar) + cofs(4) .* (v - self.vStar) >= 0) .* cofs(3);
        end

        function val = get.coefficients(self)
            val(1) = self.B .* self.windSpeed.^self.m;
            val(2) = self.view_factor .* self.emissivity .* self.stefan_boltzmann;
            val(3) = 31.9 .* self.a_SW;
            val(4) = 74.4 .* self.a_SW;
            val(5) = 5.6;
        end
    end
end