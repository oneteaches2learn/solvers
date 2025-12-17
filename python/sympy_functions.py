# newmod.py
import sympy as sp
from sympy.printing import octave_code
from sympy.parsing.sympy_parser import parse_expr

# Map MATLAB-style function names to SymPy callables
MATLAB_FUNCS = {
    # Basic trig
    'sin' :  sp.sin,
    'cos' :  sp.cos,
    'tan' :  sp.tan,
    'asin':  sp.asin,
    'acos':  sp.acos,
    'atan':  sp.atan,

    # Hyperbolic
    'sinh':  sp.sinh,
    'cosh':  sp.cosh,
    'tanh':  sp.tanh,
    'asinh': sp.asinh,
    'acosh': sp.acosh,
    'atanh': sp.atanh,

    # Exponential / log
    'exp' :  sp.exp,
    'log' :  sp.log,  # natural log

    # MATLAB log10/log2
    'log10': lambda x: sp.log(x, 10),
    'log2' : lambda x: sp.log(x, 2),

    # Roots
    'sqrt':  sp.sqrt,

    # “Scalar non-smooth” stuff
    'abs' :  sp.Abs,   # |x|
    'sign':  sp.sign,  # sign(x)

    # Min/max
    'max' :  sp.Max,
    'min' :  sp.Min,

    # Rounding
    'floor': sp.floor,
    'ceil' : sp.ceiling,  # SymPy uses "ceiling"

    # Special functions commonly used in applied math
    'erf' :  sp.erf,
    'erfc':  sp.erfc,
    'heaviside': sp.Heaviside,

    # 2-arg atan2(y,x); MATLAB uses atan2(Y,X) with same ordering
    'atan2': lambda y, x: sp.atan2(y, x),
}

# Map MATLAB-style constants to SymPy equivalents
MATLAB_CONST = {
    'pi':  sp.pi,
    'Pi':  sp.pi,
    'inf': sp.oo,
    'Inf': sp.oo,
    'INF': sp.oo,
}


def sympy_deriv(func_str,var_name="x"):

    # convert MATLAB expression to SymPy
    func_str, var_names = functionHandle2str(func_str);
    expr = matlab2sympy(func_str,var_name)

    # symbolic differentiation
    var = sp.symbols(var_name)
    df = sp.diff(expr, var)

    # convert to Octave/MATLAB code
    df_mat = octave_code(df)

    return df_mat


def matlab2sympy(func_str,var_name="x"):
    
    # replace MATLAB-style operators with Python/SymPy equivalents
    s = func_str
    s = s.replace('.^', '**')
    s = s.replace('^', '**') 
    s = s.replace('./', '/')
    s = s.replace('.*', '*')

    # build dictionary with variables, constants, and functions
    local_dict = {name: sp.symbols(name) for name in var_name}
    local_dict.update(MATLAB_CONST)
    local_dict.update(MATLAB_FUNCS)

    # parse expression
    expr = parse_expr(s, local_dict=local_dict, evaluate=True)

    return expr


def functionHandle2str(handle_str):

    s = handle_str.strip()

    # Case 1: function handle starting with '@'
    if s.startswith('@'):

        # Strip '@' and any following whitespace
        rest = s[1:].lstrip()

        # Two possible MATLAB syntaxes:
        #   @(x) ...
        #   @x ...
        if rest.startswith('('):

            # Find the closing ')'
            close_paren = rest.find(')')
            if close_paren == -1:
                raise ValueError(f"Malformed function handle (no closing ')' in {handle_str!r})")

            arglist_str = rest[1:close_paren]
            body_str    = rest[close_paren + 1:].lstrip()
            arg_names   = [a.strip() for a in arglist_str.split(',') if a.strip()]

        else:
            # Shorthand: @x x.^2  (single variable, no parentheses)
            i = 0
            while i < len(rest) and (rest[i].isalnum() or rest[i] == '_'):
                i += 1
            arg_name = rest[:i]
            body_str = rest[i:].lstrip()
            arg_names = [arg_name]

    # Case 2: just an expression, no '@'
    else:
        body_str  = s
        arg_names = list(default_symbols)    

    return body_str, arg_names
