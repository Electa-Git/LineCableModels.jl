# Function: defining functions

using .GetDPFrontend: add_raw_code, comment, make_args

"""
	Function

Defining functions.
"""
mutable struct Function <: AbstractGetDPObject
	name::String
	content::String
	comment::Union{String, Nothing}
	indent::String

	function Function()
		new("Function", "", nothing, " "^4)
	end
end

"""
	add!(func::Function, id, expression; comment=nothing)

Add a simple function to the Function object.
"""
function add!(func::Function, id, expression; comment = nothing)
	c = "$(id)() = $(expression);"
	func.content *= c

	if comment !== nothing
		add_comment!(func, comment; newline = false)
	end

	func.content *= "\n"
	return id
end

"""
	add!(func::Function, id; expression, arguments=String[], region=String[], comment=nothing)

Add a function to the Function object with optional region and arguments.
"""
function add!(
	func::Function,
	id;
	expression,
	arguments = String[],
	region = String[],
	comment = nothing,
)
	# Format the identifier
	if !isempty(region)
		region_str = make_args(region)
		id_str = "$(id)[Region[$(region_str)]]"
	else
		id_str = "$(id)[]"
	end

	# Format the expression with arguments (if any)
	expr = expression
	if !isempty(arguments)
		# Ensure expression contains placeholders for arguments
		for arg in arguments
			if !occursin(arg, expr)
				Base.error("Argument $arg not found in expression: $expr")
			end
		end
	end

	c = "$(id_str) = $(expr);"
	func.content *= c

	if comment !== nothing
		add_comment!(func, comment; newline = false)
	end

	func.content *= "\n"
	return id
end

"""
	add_list!(func::Function, id, expression_list; comment=nothing)

Add a list of functions to the Function object.
"""
function add_list!(func::Function, id, expression_list; comment = nothing)
	c = "$(id)[] = {"

	for (i, expr) in enumerate(expression_list)
		if i > 1
			c *= ", "
		end
		c *= "$(expr)"
	end

	c *= "};"
	func.content *= c

	if comment !== nothing
		add_comment!(func, comment; newline = false)
	end

	func.content *= "\n"
	return id
end

"""
	add_file!(func::Function, id, filename; comment=nothing)

Add a function from a file to the Function object.
"""
function add_file!(func::Function, id, filename; comment = nothing)
	c = "$(id)[] = Analytic[File[\"$(filename)\"]];"
	func.content *= c

	if comment !== nothing
		add_comment!(func, comment; newline = false)
	end

	func.content *= "\n"
	return id
end
"""
	add_constant!(func::Function, variable, value; comment=nothing)

Add a global constant to the Function object.
"""
function add_constant!(func::Function, variable, value; comment = nothing)
	c = "$(variable) = $(value);"
	func.content *= c

	if comment !== nothing
		add_comment!(func, comment; newline = false)
	end

	func.content *= "\n"
end

"""
	add_analytic!(func::Function, id, expression; comment=nothing)

Add an analytic function to the Function object.
"""
function add_analytic!(func::Function, id, expression; comment = nothing)
	c = "$(id)[] = Analytic[$(expression)];"
	func.content *= c

	if comment !== nothing
		add_comment!(func, comment; newline = false)
	end

	func.content *= "\n"
	return id
end

"""
	add_piecewise!(func::Function, id, x, y; comment=nothing)

Add a piecewise function to the Function object.
"""
function add_piecewise!(func::Function, id, x, y; comment = nothing)
	if length(x) != length(y)
		Base.error("x and y must have the same length")
	end

	c = "$(id)[] = InterpolationLinear[$x, $y];"
	func.content *= c

	if comment !== nothing
		add_comment!(func, comment; newline = false)
	end

	func.content *= "\n"
	return id
end

"""
	add_akima!(func::Function, id, x, y; comment=nothing)

Add an Akima interpolation function to the Function object.
"""
function add_akima!(func::Function, id, x, y; comment = nothing)
	if length(x) != length(y)
		Base.error("x and y must have the same length")
	end

	c = "$(id)[] = InterpolationAkima[$x, $y];"
	func.content *= c

	if comment !== nothing
		add_comment!(func, comment; newline = false)
	end

	func.content *= "\n"
	return id
end

"""
	add_raw_code!(func::Function, raw_code, newline=true)

Add raw code to the Function object.
"""
function add_raw_code!(func::Function, raw_code, newline = true)
	func.content = add_raw_code(func.content, raw_code, newline)
end

"""
	add_comment!(func::Function, comment_text, newline=true)

Add a comment to the Function object.
"""
function add_comment!(func::Function, comment_text, newline = true)
	add_raw_code!(func, comment(comment_text, newline = false), newline)
end
"""
	add_space!(function::Function, num_spaces::Int=1)

Add a specified number of empty lines to the function object's content for spacing in the output.
"""
function add_space!(func::Function, num_spaces::Int = 1)
	func.content *= "\n"^num_spaces
	return nothing
end
"""
	code(func::Function)

Generate GetDP code for a Function object.
"""
function code(func::Function)
	code_lines = String[]
	push!(code_lines, "\nFunction{")

	# Add the content directly since it's already formatted
	if !isempty(func.content)
		for line in split(func.content, '\n')
			if !isempty(line)
				push!(code_lines, "  " * line)
			else
				push!(code_lines, "")
			end
		end
	end

	push!(code_lines, "}")

	# Add comment if present
	if func.comment !== nothing
		return comment(func.comment) * "\n" * join(code_lines, "\n")
	else
		return join(code_lines, "\n")
	end
end

"""
	Exp(expression)

Generate the GetDP string representation for the exponential function `Exp[...]`.
Example: `Exp("x")` returns `"Exp[x]"`.
"""
function Exp(expression)
	return "Exp[$(expression)]"
end

"""
	Log(expression)

Generate the GetDP string representation for the natural logarithm function `Log[...]`.
Example: `Log("x")` returns `"Log[x]"`.
"""
function Log(expression)
	return "Log[$(expression)]"
end

"""
	Log10(expression)

Generate the GetDP string representation for the base-10 logarithm function `Log10[...]`.
Example: `Log10("x")` returns `"Log10[x]"`.
"""
function Log10(expression)
	return "Log10[$(expression)]"
end

"""
	Sqrt(expression)

Generate the GetDP string representation for the square root function `Sqrt[...]`.
Example: `Sqrt("x")` returns `"Sqrt[x]"`.
"""
function Sqrt(expression)
	return "Sqrt[$(expression)]"
end
"""
	Sin(expression)

Generate the GetDP string representation for the sine function `Sin[...]`.
"""
function Sin(expression)
	return "Sin[$(expression)]"
end

"""
	Asin(expression)

Generate the GetDP string representation for the arc sine function `Asin[...]`.
Domain: expression in [-1, 1]. Range: [-Pi/2, Pi/2]. (Real valued only).
"""
function Asin(expression)
	return "Asin[$(expression)]"
end

"""
	Cos(expression)

Generate the GetDP string representation for the cosine function `Cos[...]`.
"""
function Cos(expression)
	return "Cos[$(expression)]"
end

"""
	Acos(expression)

Generate the GetDP string representation for the arc cosine function `Acos[...]`.
Domain: expression in [-1, 1]. Range: [0, Pi]. (Real valued only).
"""
function Acos(expression)
	return "Acos[$(expression)]"
end

"""
	Tan(expression)

Generate the GetDP string representation for the tangent function `Tan[...]`.
"""
function Tan(expression)
	return "Tan[$(expression)]"
end

"""
	Atan(expression)

Generate the GetDP string representation for the arc tangent function `Atan[...]`.
Range: [-Pi/2, Pi/2]. (Real valued only).
"""
function Atan(expression)
	return "Atan[$(expression)]"
end

"""
	Atan2(y, x)

Generate the GetDP string representation for the arc tangent function `Atan2[y, x]`.
Computes atan(y/x). Range: [-Pi, Pi]. (Real valued only).
"""
function Atan2(y, x)
	return "Atan2[$(y),$(x)]"
end

# Hyperbolic Functions
"""
	Sinh(expression)

Generate the GetDP string representation for the hyperbolic sine function `Sinh[...]`.
"""
function Sinh(expression)
	return "Sinh[$(expression)]"
end

"""
	Cosh(expression)

Generate the GetDP string representation for the hyperbolic cosine function `Cosh[...]`.
"""
function Cosh(expression)
	return "Cosh[$(expression)]"
end

"""
	Tanh(expression)

Generate the GetDP string representation for the hyperbolic tangent function `Tanh[...]`.
(Real valued only).
"""
function Tanh(expression)
	return "Tanh[$(expression)]"
end

"""
	TanhC2(expression)

Generate the GetDP string representation for the complex hyperbolic tangent function `TanhC2[...]`.
"""
function TanhC2(expression)
	return "TanhC2[$(expression)]"
end

# Basic Math/Rounding Functions
"""
	Fabs(expression)

Generate the GetDP string representation for the absolute value function `Fabs[...]`.
(Real valued only).
"""
function Fabs(expression)
	return "Fabs[$(expression)]"
end

"""
	Abs(expression)

Generate the GetDP string representation for the absolute value/modulus function `Abs[...]`.
(Works for complex numbers).
"""
function Abs(expression)
	return "Abs[$(expression)]"
end

"""
	Floor(expression)

Generate the GetDP string representation for the floor function `Floor[...]`.
Rounds downwards. (Real valued only).
"""
function Floor(expression)
	return "Floor[$(expression)]"
end

"""
	Ceil(expression)

Generate the GetDP string representation for the ceiling function `Ceil[...]`.
Rounds upwards. (Real valued only).
"""
function Ceil(expression)
	return "Ceil[$(expression)]"
end

"""
	Fmod(x, y)

Generate the GetDP string representation for the floating-point remainder function `Fmod[x, y]`.
Remainder of x/y with sign of x. (Real valued only).
"""
function Fmod(x, y)
	return "Fmod[$(x),$(y)]"
end

"""
	Min(a, b)

Generate the GetDP string representation for the minimum function `Min[a, b]`.
(Scalar, real valued only).
"""
function Min(a, b)
	return "Min[$(a),$(b)]"
end

"""
	Max(a, b)

Generate the GetDP string representation for the maximum function `Max[a, b]`.
(Scalar, real valued only).
"""
function Max(a, b)
	return "Max[$(a),$(b)]"
end

"""
	Sign(expression)

Generate the GetDP string representation for the sign function `Sign[...]`.
Returns -1 for expression < 0, 1 otherwise. (Real valued only).
"""
function Sign(expression)
	return "Sign[$(expression)]"
end

# Bessel Functions (assuming GetDP syntax Jn[order, value])
"""
	Jn(order, value)

Generate the GetDP string representation for the Bessel function of the first kind `Jn[order, value]`.
(Real valued only).
"""
function Jn(order, value)
	return "Jn[$(order),$(value)]"
end

"""
	dJn(order, value)

Generate the GetDP string representation for the derivative of the Bessel function of the first kind `dJn[order, value]`.
(Real valued only).
"""
function dJn(order, value)
	return "dJn[$(order),$(value)]"
end

"""
	Yn(order, value)

Generate the GetDP string representation for the Bessel function of the second kind `Yn[order, value]`.
(Real valued only).
"""
function Yn(order, value)
	return "Yn[$(order),$(value)]"
end

"""
	dYn(order, value)

Generate the GetDP string representation for the derivative of the Bessel function of the second kind `dYn[order, value]`.
(Real valued only).
"""
function dYn(order, value)
	return "dYn[$(order),$(value)]"
end


# Vector/Tensor Functions
"""
	Cross(vector1, vector2)

Generate the GetDP string representation for the cross product `Cross[vector1, vector2]`.
Arguments must be vectors.
"""
function Cross(vector1, vector2)
	return "Cross[$(vector1),$(vector2)]"
end

"""
	Hypot(a, b)

Generate the GetDP string representation for the hypotenuse function `Hypot[a, b]`.
Computes Sqrt(a^2 + b^2).
"""
function Hypot(a, b)
	return "Hypot[$(a),$(b)]"
end

"""
	Norm(expression)

Generate the GetDP string representation for the norm function `Norm[...]`.
Absolute value for scalar, Euclidean norm for vector.
"""
function Norm(expression)
	return "Norm[$(expression)]"
end

"""
	SquNorm(expression)

Generate the GetDP string representation for the square norm function `SquNorm[...]`.
Equivalent to Norm[expression]^2.
"""
function SquNorm(expression)
	return "SquNorm[$(expression)]"
end

"""
	Unit(expression)

Generate the GetDP string representation for the unit vector function `Unit[...]`.
Computes expression / Norm[expression]. Returns 0 if norm is near zero.
"""
function Unit(expression)
	return "Unit[$(expression)]"
end

"""
	Transpose(expression)

Generate the GetDP string representation for the transpose function `Transpose[...]`.
Expression must be a tensor.
"""
function Transpose(expression)
	return "Transpose[$(expression)]"
end

"""
	Inv(expression)

Generate the GetDP string representation for the inverse tensor function `Inv[...]`.
Expression must be a tensor.
"""
function Inv(expression)
	return "Inv[$(expression)]"
end

"""
	Det(expression)

Generate the GetDP string representation for the tensor determinant function `Det[...]`.
Expression must be a tensor.
"""
function Det(expression)
	return "Det[$(expression)]"
end

"""
	Rotate(object, rot_x, rot_y, rot_z)

Generate the GetDP string representation for the rotation function `Rotate[object, rx, ry, rz]`.
Rotates a vector or tensor `object` by angles `rot_x`, `rot_y`, `rot_z` (radians) around axes x, y, z.
"""
function Rotate(object_to_rotate, rot_x, rot_y, rot_z)
	return "Rotate[$(object_to_rotate),$(rot_x),$(rot_y),$(rot_z)]"
end

"""
	TTrace(tensor)

Generate the GetDP string representation for the tensor trace function `TTrace[...]`.
Expression must be a tensor.
"""
function TTrace(tensor)
	return "TTrace[$(tensor)]"
end

# Special/Time Functions
"""
	Cos_wt_p(omega, phase)

Generate the GetDP string representation for the time function `Cos_wt_p[]{omega, phase}`.
Real: Cos[omega*Time + phase]. Complex: Complex[Cos[phase], Sin[phase]].
"""
function Cos_wt_p(omega, phase)
	# Note the empty [] and double {{}}
	return "Cos_wt_p[]{{$(omega),$(phase)}}"
end

"""
	Sin_wt_p(omega, phase)

Generate the GetDP string representation for the time function `Sin_wt_p[]{omega, phase}`.
Real: Sin[omega*Time + phase]. Complex: Complex[Sin[phase], -Cos[phase]].
"""
function Sin_wt_p(omega, phase)
	# Note the empty [] and double {{}}
	return "Sin_wt_p[]{{$(omega),$(phase)}}"
end

"""
	Period(expression, period_const)

Generate the GetDP string representation for the periodic function `Period[expr]{period_const}`.
Result is always in [0, period_const[.
"""
function _period_expression(expression, period_const)
	# Note the {} after []
	return "Period[$(expression)]{$(period_const)}"
end
# Green Functions
"""
	Laplace(dim)

Generate the GetDP string representation for the Laplace Green function `Laplace[]{dim}`.
Result depends on `dim`: r/2 (1D), (1/(2*Pi))*ln(1/r) (2D), 1/(4*Pi*r) (3D).
"""
function Laplace(dim)
	return "Laplace[]{{{$(dim)}}}" # Note: empty [], triple {{{}}}
end

"""
	GradLaplace(dim)

Generate the GetDP string representation for the gradient of the Laplace Green function `GradLaplace[]{dim}`.
Gradient is relative to the destination point (X, Y, Z).
"""
function GradLaplace(dim)
	return "GradLaplace[]{{{$(dim)}}}" # Note: empty [], triple {{{}}}
end

"""
	Helmholtz(dim, k0)

Generate the GetDP string representation for the Helmholtz Green function `Helmholtz[]{dim, k0}`.
Typically exp(j*k0*r)/(4*Pi*r) (3D). `k0` is the wave number.
"""
function Helmholtz(dim, k0)
	return "Helmholtz[]{{{$(dim),$(k0)}}}" # Note: empty [], triple {{{}}}
end

"""
	GradHelmholtz(dim, k0)

Generate the GetDP string representation for the gradient of the Helmholtz Green function `GradHelmholtz[]{dim, k0}`.
Gradient is relative to the destination point (X, Y, Z).
"""
function GradHelmholtz(dim, k0)
	return "GradHelmholtz[]{{{$(dim),$(k0)}}}" # Note: empty [], triple {{{}}}
end

# Type Manipulation Functions
"""
	Complex(re1, im1, re2, im2, ...)

Generate the GetDP string representation for creating a complex value `Complex[re1, im1, ...]`.
Takes an even number of real-valued expressions.
"""
function _complex_expression(args...) # Varargs
	arg_string = join(args, ',')
	return "Complex[$(arg_string)]"
end

"""
	Re(complex_expr)

Generate the GetDP string representation for the real part function `Re[...]`.
"""
function Re(complex_expr)
	return "Re[$(complex_expr)]"
end

"""
	Im(complex_expr)

Generate the GetDP string representation for the imaginary part function `Im[...]`.
"""
function Im(complex_expr)
	return "Im[$(complex_expr)]"
end

"""
	Conj(complex_expr)

Generate the GetDP string representation for the complex conjugate function `Conj[...]`.
"""
function Conj(complex_expr)
	return "Conj[$(complex_expr)]"
end

"""
	Cart2Pol(complex_expr)

Generate the GetDP string representation for Cartesian to Polar conversion `Cart2Pol[...]`.
Input: Complex[real, imag]. Output: Complex[amplitude, phase].
"""
function Cart2Pol(complex_expr)
	return "Cart2Pol[$(complex_expr)]"
end


# Vector/Tensor Creation
"""
	Vector(s0, s1, s2)

Generate the GetDP string representation for creating a vector `Vector[s0, s1, s2]`.
"""
function _vector_expression(s0, s1, s2)
	return "Vector[$(s0),$(s1),$(s2)]"
end

"""
	Tensor(s00, s01, s02, s10, s11, s12, s20, s21, s22)

Generate the GetDP string representation for creating a tensor from 9 scalars (row-major) `Tensor[...]`.
"""
function Tensor(s00, s01, s02, s10, s11, s12, s20, s21, s22)
	return "Tensor[$(s00),$(s01),$(s02),$(s10),$(s11),$(s12),$(s20),$(s21),$(s22)]"
end

"""
	TensorV(v0, v1, v2)

Generate the GetDP string representation for creating a tensor from 3 row vectors `TensorV[...]`.
"""
function TensorV(v0, v1, v2)
	return "TensorV[$(v0),$(v1),$(v2)]"
end

"""
	TensorSym(s00, s01, s02, s11, s12, s22)

Generate the GetDP string representation for creating a symmetric tensor from 6 scalars `TensorSym[...]`.
T = [s00 s01 s02; s01 s11 s12; s02 s12 s22]
"""
function TensorSym(s00, s01, s02, s11, s12, s22)
	return "TensorSym[$(s00),$(s01),$(s02),$(s11),$(s12),$(s22)]"
end

"""
	TensorDiag(s00, s11, s22)

Generate the GetDP string representation for creating a diagonal tensor from 3 scalars `TensorDiag[...]`.
"""
function TensorDiag(s00, s11, s22)
	return "TensorDiag[$(s00),$(s11),$(s22)]"
end

"""
	SquDyadicProduct(vector)

Generate the GetDP string representation for the square dyadic product `SquDyadicProduct[...]`.
Computes vector * Transpose[vector].
"""
function SquDyadicProduct(vector)
	return "SquDyadicProduct[$(vector)]"
end

# Component Extraction
"""
	CompX(vector)

Generate the GetDP string representation for getting the X component `CompX[...]`.
"""
function CompX(vector)
	return "CompX[$(vector)]"
end

"""
	CompY(vector)

Generate the GetDP string representation for getting the Y component `CompY[...]`.
"""
function CompY(vector)
	return "CompY[$(vector)]"
end

"""
	CompZ(vector)

Generate the GetDP string representation for getting the Z component `CompZ[...]`.
"""
function CompZ(vector)
	return "CompZ[$(vector)]"
end

"""
	CompXX(tensor)

Generate the GetDP string representation for getting the XX component `CompXX[...]`.
"""
function CompXX(tensor)
	return "CompXX[$(tensor)]"
end

"""
	CompYY(tensor)

Generate the GetDP string representation for getting the YY component `CompYY[...]`.
"""
function CompYY(tensor)
	return "CompYY[$(tensor)]"
end

"""
	CompZZ(tensor)

Generate the GetDP string representation for getting the ZZ component `CompZZ[...]`.
"""
function CompZZ(tensor)
	return "CompZZ[$(tensor)]"
end

"""
	CompXY(tensor)

Generate the GetDP string representation for getting the XY component `CompXY[...]`.
"""
function CompXY(tensor)
	return "CompXY[$(tensor)]"
end

"""
	CompYX(tensor)

Generate the GetDP string representation for getting the YX component `CompYX[...]`.
"""
function CompYX(tensor)
	return "CompYX[$(tensor)]"
end

"""
	CompXZ(tensor)

Generate the GetDP string representation for getting the XZ component `CompXZ[...]`.
"""
function CompXZ(tensor)
	return "CompXZ[$(tensor)]"
end

"""
	CompZX(tensor)

Generate the GetDP string representation for getting the ZX component `CompZX[...]`.
"""
function CompZX(tensor)
	return "CompZX[$(tensor)]"
end

"""
	CompYZ(tensor)

Generate the GetDP string representation for getting the YZ component `CompYZ[...]`.
"""
function CompYZ(tensor)
	return "CompYZ[$(tensor)]"
end

"""
	CompZY(tensor)

Generate the GetDP string representation for getting the ZY component `CompZY[...]`.
"""
function CompZY(tensor)
	return "CompZY[$(tensor)]"
end


# Coordinate Transformations
"""
	Cart2Sph(vector)

Generate the GetDP string representation for the Cartesian to Spherical transformation tensor `Cart2Sph[...]`.
"""
function Cart2Sph(vector)
	return "Cart2Sph[$(vector)]"
end

"""
	Cart2Cyl(vector)

Generate the GetDP string representation for the Cartesian to Cylindrical transformation tensor `Cart2Cyl[...]`.
"""
function Cart2Cyl(vector)
	return "Cart2Cyl[$(vector)]"
end

# Unit Vectors
"""
	UnitVectorX()

Generate the GetDP string representation for the unit vector in X: `UnitVectorX[]`.
"""
function UnitVectorX()
	return "UnitVectorX[]"
end

"""
	UnitVectorY()

Generate the GetDP string representation for the unit vector in Y: `UnitVectorY[]`.
"""
function UnitVectorY()
	return "UnitVectorY[]"
end

"""
	UnitVectorZ()

Generate the GetDP string representation for the unit vector in Z: `UnitVectorZ[]`.
"""
function UnitVectorZ()
	return "UnitVectorZ[]"
end


# Coordinate Functions
"""
	X()

Generate the GetDP string representation for the X coordinate: `X[]`.
"""
function X()
	return "X[]"
end

"""
	Y()

Generate the GetDP string representation for the Y coordinate: `Y[]`.
"""
function Y()
	return "Y[]"
end

"""
	Z()

Generate the GetDP string representation for the Z coordinate: `Z[]`.
"""
function Z()
	return "Z[]"
end

"""
	XYZ()

Generate the GetDP string representation for the coordinate vector: `XYZ[]`.
"""
function XYZ()
	return "XYZ[]"
end
# Miscellaneous Functions
"""
	Printf(expression)

Generate the GetDP string representation for printing a value during evaluation: `Printf[expression]`.
"""
function Printf(expression)
	return "Printf[$(expression)]"
end

"""
	Rand(max_val)

Generate the GetDP string representation for a pseudo-random number in [0, max_val]: `Rand[max_val]`.
"""
function Rand(max_val)
	return "Rand[$(max_val)]"
end

"""
	Normal()

Generate the GetDP string representation for the element's normal vector: `Normal[]`.
"""
function Normal()
	return "Normal[]"
end

"""
	NormalSource()

Generate the GetDP string representation for the source element's normal vector: `NormalSource[]`.
(Valid in Integral quantity).
"""
function NormalSource()
	return "NormalSource[]"
end

"""
	Tangent()

Generate the GetDP string representation for the element's tangent vector: `Tangent[]`.
(Valid for line elements).
"""
function Tangent()
	return "Tangent[]"
end

"""
	TangentSource()

Generate the GetDP string representation for the source element's tangent vector: `TangentSource[]`.
(Valid in Integral quantity, line elements).
"""
function TangentSource()
	return "TangentSource[]"
end

"""
	ElementVol()

Generate the GetDP string representation for the element's volume (or area/length): `ElementVol[]`.
"""
function ElementVol()
	return "ElementVol[]"
end

"""
	SurfaceArea(list_expr::String="")

Generate the GetDP string representation for surface area calculation: `SurfaceArea[]{list}`.
`list_expr` is a comma-separated string of physical surface tags, or empty for the current surface.
"""
function SurfaceArea(list_expr::String = "")
	return "SurfaceArea[]{{{$(list_expr)}}}" # Note: empty [], triple {{{}}}
end

"""
	GetVolume()

Generate the GetDP string representation for the volume of the current physical group: `GetVolume[]`.
"""
function GetVolume()
	return "GetVolume[]"
end

"""
	CompElementNum()

Generate the GetDP string representation to compare current and source element tags: `CompElementNum[]`.
Returns 0 if identical.
"""
function CompElementNum()
	return "CompElementNum[]"
end

"""
	GetNumElements(list_expr::String="")

Generate the GetDP string representation for counting elements: `GetNumElements[]{list}`.
`list_expr` is a comma-separated string of physical region tags, or empty for the current region.
"""
function GetNumElements(list_expr::String = "")
	return "GetNumElements[]{{{$(list_expr)}}}" # Note: empty [], triple {{{}}}
end

"""
	ElementNum()

Generate the GetDP string representation for the current element's tag: `ElementNum[]`.
"""
function ElementNum()
	return "ElementNum[]"
end

"""
	QuadraturePointIndex()

Generate the GetDP string representation for the current quadrature point index: `QuadraturePointIndex[]`.
"""
function QuadraturePointIndex()
	return "QuadraturePointIndex[]"
end

"""
	AtIndex(index_expr, list_expr::String)

Generate the GetDP string representation for accessing list element by index: `AtIndex[index]{list}`.
`list_expr` is a comma-separated string list. Index is 0-based(? check GetDP docs).
"""
function AtIndex(index_expr, list_expr::String)
	return "AtIndex[$(index_expr)]{$(list_expr)}" # Note: []{}
end

# Interpolation Functions
"""
	InterpolationLinear(x_expr, table_list_expr::String)

Generate the GetDP string representation for linear interpolation: `InterpolationLinear[x]{table}`.
`table_list_expr` is a comma-separated list of x,y pairs: "x1,y1,x2,y2,...".
"""
function InterpolationLinear(x_expr, table_list_expr::String)
	return "InterpolationLinear[$(x_expr)]{$(table_list_expr)}" # Note: []{}
end

"""
	dInterpolationLinear(x_expr, table_list_expr::String)

Generate the GetDP string representation for the derivative of linear interpolation: `dInterpolationLinear[x]{table}`.
"""
function dInterpolationLinear(x_expr, table_list_expr::String)
	return "dInterpolationLinear[$(x_expr)]{$(table_list_expr)}" # Note: []{}
end

"""
	InterpolationBilinear(x_expr, y_expr, table_list_expr::String)

Generate the GetDP string representation for bilinear interpolation: `InterpolationBilinear[x, y]{table}`.
Table format needs checking in GetDP docs.
"""
function InterpolationBilinear(x_expr, y_expr, table_list_expr::String)
	return "InterpolationBilinear[$(x_expr), $(y_expr)]{$(table_list_expr)}" # Note: []{}
end

"""
	dInterpolationBilinear(x_expr, y_expr, table_list_expr::String)

Generate the GetDP string representation for the derivative of bilinear interpolation: `dInterpolationBilinear[x, y]{table}`.
Result is a vector.
"""
function dInterpolationBilinear(x_expr, y_expr, table_list_expr::String)
	return "dInterpolationBilinear[$(x_expr), $(y_expr)]{$(table_list_expr)}" # Note: []{}
end

"""
	InterpolationAkima(x_expr, table_list_expr::String)

Generate the GetDP string representation for Akima interpolation: `InterpolationAkima[x]{table}`.
`table_list_expr` is a comma-separated list of x,y pairs: "x1,y1,x2,y2,...".
"""
function InterpolationAkima(x_expr, table_list_expr::String)
	return "InterpolationAkima[$(x_expr)]{$(table_list_expr)}" # Note: []{}
end

"""
	dInterpolationAkima(x_expr, table_list_expr::String)

Generate the GetDP string representation for the derivative of Akima interpolation: `dInterpolationAkima[x]{table}`.
"""
function dInterpolationAkima(x_expr, table_list_expr::String)
	return "dInterpolationAkima[$(x_expr)]{$(table_list_expr)}" # Note: []{}
end

"""
	Order(quantity_name)

Generate the GetDP string representation for getting interpolation order: `Order[quantity]`.
"""
function Order(quantity_name)
	return "Order[$(quantity_name)]"
end


# Field Evaluation Functions
"""
	Field(eval_point_expr)

Generate the GetDP string representation for evaluating the last Gmsh field: `Field[eval_point]`.
Typically `eval_point_expr` is `XYZ[]`.
"""
function Field(eval_point_expr)
	return "Field[$(eval_point_expr)]"
end

"""
	Field(eval_point_expr, tags_list_expr::String)

Generate the GetDP string representation for evaluating and summing specific Gmsh fields: `Field[eval_point]{tags_list}`.
`tags_list_expr` is a comma-separated list of field tags.
"""
function Field(eval_point_expr, tags_list_expr::String)
	return "Field[$(eval_point_expr)]{$(tags_list_expr)}" # Note: []{}
end

# Common helper for typed field functions
function _TypedFieldHelper(
	func_name,
	expression,
	expression_cst_list,
	timestep,
	elmt_interp,
)
	interp_flag = Int(elmt_interp) # Convert Bool to 0 or 1
	return "$(func_name)[$(expression), $(timestep), $(interp_flag)]{$(expression_cst_list)}"
end

"""
	ScalarField(expression, expression_cst_list; timestep=0, elmt_interp=true)

Generate GetDP string for evaluating scalar fields: `ScalarField[expr, ts, interp]{list}`.
"""
function ScalarField(
	expression,
	expression_cst_list::String;
	timestep = 0,
	elmt_interp::Bool = true,
)
	return _TypedFieldHelper(
		"ScalarField",
		expression,
		expression_cst_list,
		timestep,
		elmt_interp,
	)
end

"""
	VectorField(expression, expression_cst_list; timestep=0, elmt_interp=true)

Generate GetDP string for evaluating vector fields: `VectorField[expr, ts, interp]{list}`.
"""
function VectorField(
	expression,
	expression_cst_list::String;
	timestep = 0,
	elmt_interp::Bool = true,
)
	return _TypedFieldHelper(
		"VectorField",
		expression,
		expression_cst_list,
		timestep,
		elmt_interp,
	)
end

"""
	TensorField(expression, expression_cst_list; timestep=0, elmt_interp=true)

Generate GetDP string for evaluating tensor fields: `TensorField[expr, ts, interp]{list}`.
"""
function TensorField(
	expression,
	expression_cst_list::String;
	timestep = 0,
	elmt_interp::Bool = true,
)
	return _TypedFieldHelper(
		"TensorField",
		expression,
		expression_cst_list,
		timestep,
		elmt_interp,
	)
end

"""
	ComplexScalarField(expression, expression_cst_list; timestep=0, elmt_interp=true)

Generate GetDP string for evaluating complex scalar fields: `ComplexScalarField[expr, ts, interp]{list}`.
"""
function ComplexScalarField(
	expression,
	expression_cst_list::String;
	timestep = 0,
	elmt_interp::Bool = true,
)
	return _TypedFieldHelper(
		"ComplexScalarField",
		expression,
		expression_cst_list,
		timestep,
		elmt_interp,
	)
end

"""
	ComplexVectorField(expression, expression_cst_list; timestep=0, elmt_interp=true)

Generate GetDP string for evaluating complex vector fields: `ComplexVectorField[expr, ts, interp]{list}`.
"""
function ComplexVectorField(
	expression,
	expression_cst_list::String;
	timestep = 0,
	elmt_interp::Bool = true,
)
	return _TypedFieldHelper(
		"ComplexVectorField",
		expression,
		expression_cst_list,
		timestep,
		elmt_interp,
	)
end

"""
	ComplexTensorField(expression, expression_cst_list; timestep=0, elmt_interp=true)

Generate GetDP string for evaluating complex tensor fields: `ComplexTensorField[expr, ts, interp]{list}`.
"""
function ComplexTensorField(
	expression,
	expression_cst_list::String;
	timestep = 0,
	elmt_interp::Bool = true,
)
	return _TypedFieldHelper(
		"ComplexTensorField",
		expression,
		expression_cst_list,
		timestep,
		elmt_interp,
	)
end


# Runtime Information/Control Functions
"""
	GetCpuTime()

Generate GetDP string for getting CPU time: `GetCpuTime[]`.
"""
function GetCpuTime()
	return "GetCpuTime[]"
end

"""
	GetWallClockTime()

Generate GetDP string for getting wall clock time: `GetWallClockTime[]`.
"""
function GetWallClockTime()
	return "GetWallClockTime[]"
end

"""
	GetMemory()

Generate GetDP string for getting memory usage (MB): `GetMemory[]`.
"""
function GetMemory()
	return "GetMemory[]"
end

"""
	SetNumberRunTime(value_expr, name::String)

Generate GetDP string for setting a ONELAB variable at runtime: `SetNumberRunTime[value]{"name"}`.
"""
function SetNumberRunTime(value_expr, name::String)
	# Note the escaped quotes inside {}
	return "SetNumberRunTime[$(value_expr)]{\"$name\"}"
end

"""
	GetNumberRunTime(name::String; default_value=nothing)

Generate GetDP string for getting a ONELAB variable at runtime: `GetNumberRunTime["name"]` or `GetNumberRunTime["name"]{default}`.
"""
function GetNumberRunTime(name::String; default_value = nothing)
	# Note the escaped quotes inside []
	if default_value === nothing
		return "GetNumberRunTime[\"$name\"]"
	else
		return "GetNumberRunTime[\"$name\"]{$(default_value)}"
	end
end

"""
	SetVariable(value_expr, variable_id::String)

Generate GetDP string for setting a runtime variable: `SetVariable[value]{variable_id}`.

"""
function SetVariable(value_expr, variable_id::String)
	return "SetVariable[$(value_expr)]{\$$(variable_id)}" # Note: \$ escapes $
end

"""
	GetVariable(variable_id::String; default_value=nothing)

Generate GetDP string for getting a runtime variable: `GetVariable[]{variable_id}` or `GetVariable[default]{variable_id}`.
"""
function GetVariable(variable_id::String; default_value = nothing)
	# Note: \$ escapes $
	if default_value === nothing
		return "GetVariable[]{\$$(variable_id)}"
	else
		return "GetVariable[$(default_value)]{\$$(variable_id)}"
	end
end


# Index/Table Functions
"""
	ValueFromIndex(list_expr::String)

Generate GetDP string for getting value from index map: `ValueFromIndex[]{list}`.
List format: "entity1, value1, entity2, value2, ...".
"""
function ValueFromIndex(list_expr::String)
	return "ValueFromIndex[]{{{$(list_expr)}}}" # Note: empty [], triple {{{}}}
end

"""
	VectorFromIndex(list_expr::String)

Generate GetDP string for getting vector from index map: `VectorFromIndex[]{list}`.
List format: "entity1, v1x, v1y, v1z, entity2, v2x, ...".
"""
function VectorFromIndex(list_expr::String)
	return "VectorFromIndex[]{{{$(list_expr)}}}" # Note: empty [], triple {{{}}}
end

"""
	ValueFromTable(default_expr, table_name::String)

Generate GetDP string for getting value from PostOperation table: `ValueFromTable[default]{"table_name"}`.
"""
function ValueFromTable(default_expr, table_name::String)
	# Note escaped quotes inside {}
	return "ValueFromTable[$(default_expr)]{\"$table_name\"}"
end
