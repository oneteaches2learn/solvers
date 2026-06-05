classdef GaussianFunction
    properties
        expr sym
        vars sym
        params struct
    end

    methods
        function self = GaussianFunction(expr, vars)
            self.expr = expr;
            self.vars = vars;
            self.params = struct();
        end

        function self = addParameter(self, name, param)
            arguments
                self 
                name string
                param
            end

            self.params.(name) = param;
        end

        function [fh, sampledParams] = sample(self)

            sampledParams = struct();

            exprSampled = self.expr;
            names = fieldnames(self.params);

            for k = 1:numel(names)

                name = names{k};
                p = self.params.(name);

                if isa(p, "GaussianParameter")
                    value = p.sample();
                elseif isnumeric(p)
                    value = p;
                else
                    error("GaussianSymbolicFunction:InvalidParameter", ...
                          "Parameter %s must be numeric or GaussianParameter.", name);
                end

                sampledParams.(name) = value;
                exprSampled = subs(exprSampled, sym(name), value);

            end

            fh = matlabFunction(exprSampled);
        end
    end
end