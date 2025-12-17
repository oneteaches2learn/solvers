classdef ManufacturedFunctions_ode < Coefficients
% ManufacturedFunctions2d_ode manufactures the source for an ODE MMS test.
%
% The simplest ODE problem is of the form
%
%			c dv/dt + s(u,v) = g
%
% In an MMS test, vTrue and the coefficients c is presumed known.
% Source g is 'manufactured' from these quantities. The ManufacturedFunctions_ode
% class manufactures and stores g. This class also manufactures ds/dv for Newton
% solvers and the initial condition vInit.
% 
% Author: Tyler Fara			Date: Nov 11, 2025
%-----------------------------------------------------------------------------%
% Notes
%  (1) ManufacturedFunctions_ode takes advantage of MATLAB's symbolic functions
%  toolkit to compute derivatives and to assemble the resulting source function
%  g. As a result, the inputs to the ManufacturedFunctions_ode constructor and
%  all of the properties of ManufacturedFunctions_ode are symfun objects.
%
%  (2) Ultimately, the manufactured source g, the coefficients c and s, as well
%  as ds/dv will be given to a solver, and solvers take function_handle objects.
%  Thus, the ManufacturedFunctions_ode class contains methods to convert each
%  symfun to its corresponding function_handle. 
%
%  (3) The output of the coefficients2FunctionHandles method is 
%  a structure, each property of which is a function_handle. 
%
%  (4) This class manufactures everything automatically using getters. If this
%  idea works, then you should go back to your previous ManufacturedFunctions2d
%  classes and refactor them to use getters as well.
% -----------------------------------------------------------------------------%


	properties
        vTrue   % symfun, true solution vTrue = vTrue(t)
        uTrue   % symfun, true solution uTrue = uTrue(x,t)
        uTrue_avg % symfun, uTrue integrated over the spatial domain
		c		% double or symfun, coefficient in c dv/dt
        c_t     % double or symfun, time derivative of c, i.e. c_t = dc/dt
		vInit	% symfun, vInit(t) = vTrue(0)
		v_t 	% symfun, v_t = d(vTrue)/dt, i.e. is the time derivative of vTrue 
		cv_t	% symfun, cv_t = d(c * vTrue)/dt, i.e
        s       % symfun, s = s(u,v)
        ds      % symfun, ds = ds/dv
        g		% symfun, manufactured source function
	end

	methods
		function self = ManufacturedFunctions_ode(c,s,vTrue,uTrue)
		% ManufacturedFunctions2d(p,k,uTrue) inputs are symfun objects
			
			% call superclass constructor
			self@Coefficients();

			% make vTrue, uTrue function of time and time and space, respectively
			x = sym('x',[1 2]); syms t;
			self.vTrue = symfun(vTrue,t);
            self.uTrue = symfun(uTrue,[x,t]);

			% store input data
			self.c = c;
            self.s = s;
             
		end

        function val = get.uTrue_avg(self)
        % presumes that spatial domain is Omega = [0,1] x [0,1]

            % computes uTrue integrated over the spatial domain
            x = sym('x',[1 2]); syms t;
            val = int(int(self.uTrue,x(1),0,1),x(2),0,1);

        end

        function val = get.vInit(self)
        % evaluates vInit at t = 0
            syms t;
            val = double(self.vTrue(0));
        end

        function val = get.ds(self)
        % manufactures ds/dv
            syms t u v;
            val = diff(self.s,v);
        end

        function val = get.v_t(self)
        % manufactures v_t, the time derivative of vTrue
            syms t;
            val = diff(self.vTrue,t);
        end

        function val = get.cv_t(self)
        % manufactures cv_t, the time derivative of c * vTrue
            syms t;
            val = diff(self.c * self.vTrue,t);
        end

        function val = get.g(self)
        % manufactures g, the right-hand side function

            % insert vTrue and uTrue_avg into s(u,v)
            syms t u v;
            temp = symfun(self.s,[u,v]);
            temp = compose(temp,self.vTrue);
            temp = compose(temp,self.uTrue_avg);

            % assemble g
            val = self.cv_t + temp;

        end

		function data = functionHandles(self)

			% output struct of function_handles
			syms t u v;

            % source function as function_handle
            data.g = matlabFunction(symfun(self.g,t));

            % initial condition as scalar
			data.vInit = self.vInit;

            % coefficients as function_handles
			data.cofs.c = matlabFunction(symfun(self.c,t));
            data.cofs.s = matlabFunction(symfun(self.s,[u,v]));
            data.cofs.ds = matlabFunction(symfun(self.ds,[u,v]));

            % temp: overwrite coefficient c as a scalar
            data.cofs.c = self.c;

		end
	end
    

    % static methods
    methods (Static)

        function data = addTimeData(data,T,dt)
        % adds time data to data struct
            data.time.T = T;
            data.time.dt = dt;
        end

        function data = addConstraintData(data,constraints)
        % adds constraint data to data struct
            data.constraints = constraints;
        end
    end

end
