% sympy_test.m

fprintf('Running sympy_test.m\n');

% reload the module to get latest changes
clear classes
mod = py.importlib.import_module('newmod');
py.importlib.reload(mod);

% call the python function
func = "sin(pi/4 * x) .* exp(3*x)";
func = "sqrt(x) + log(x + 1)";
func = @(x) sqrt(x) + log(x + 1);

func = func2str(func);
df = py.newmod.sympy_deriv(func, "x");

% convert from python string to matlab string
df = string(df);

% build function handle from string
df_mat = str2func("@(x) " + df);

% test result
df_mat
x = linspace(0, 1, 5);
y = df_mat(x)


%{
% Ensure Python environment is set up
if count(py.sys.path,pwd) == 0
    insert(py.sys.path,int32(0),pwd);
end

% Example: f(x) = x^2 * sin(x), written in Python syntax
expr_py = "x**2*sp.sin(x)";     % could also do "x**2*sympy.sin(x)" etc.

% Call Python to get MATLAB-style derivative expression
mat_expr_py = py.sympy_deriv.derivative_matlab_expr(expr_py, "x");
mat_expr    = string(mat_expr_py);   % convert to MATLAB string

% Build a function handle @(x) <mat_expr>
fh = str2func("@(x) " + mat_expr);

% Test
x = linspace(0, 1, 5);
y = fh(x);
%}
