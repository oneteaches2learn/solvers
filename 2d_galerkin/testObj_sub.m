classdef testObj_sub < testObj_super

	properties
	end

	methods
		function self = testObj_super
		end

		function val = test_method(self)
			%val = 71;
			val = test_method@testObj_super(self);
		end
	end

end
