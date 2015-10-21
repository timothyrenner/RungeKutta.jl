module RungeKutta

# Implements a Runge-Kutta solver for systems of ODEs.
# Currently only the fourth order fixed width version is implemented.

# Each of these functions defines the constants for the Runge-Kutta iteration.
k1{T<:AbstractFloat}(f::Array{Function,1}, x::Array{T,1},
    t::AbstractFloat, h::AbstractFloat) =
    map(fe -> fe(t,x), f);

k2{T<:AbstractFloat}(f::Array{Function,1}, x::Array{T,1}, t::AbstractFloat,
    h::AbstractFloat, k1val::Array{T,1}) =
    map(fe -> fe(t + h/2, x + (0.5*h).*k1val), f);

k3{T<:AbstractFloat}(f::Array{Function,1}, x::Array{T,1}, t::AbstractFloat,
    h::AbstractFloat, k2val::Array{T,1}) =
    map(fe -> fe(t + h/2, x + (0.5*h).*k2val), f);

k4{T<:AbstractFloat}(f::Array{Function,1}, x::Array{T,1}, t::AbstractFloat,
    h::AbstractFloat, k3val::Array{T,1}) =
   map(fe -> fe(t + h, x + h.*k3val), f);

# Computes the value of the next point in the Runge-Kutta iteration.
#
# Args:
#   f: The array of functions for computing x_n+1.
#   x: The current value of the points in the space.
#   t: The current value of time.
#   h: The step size.
#
# Returns: The next point in the solution.
function nextPoint{T<:AbstractFloat}(f::Array{Function,1}, x::Array{T,1},
    t::AbstractFloat, h::AbstractFloat)

    #Pre-calculate the constants.
    k1val = k1(f, x, t, h);
    k2val = k2(f, x, t, h, k1val);
    k3val = k3(f, x, t, h, k2val);
    k4val = k4(f, x, t, h, k3val);

   #This is the basic Runge-Kutta fixed width method.
   return t + h, x + (h/6).*(k1val + 2k2val + 2k3val + k4val);

end #Close nextPoint.

# Solves a system of ODEs using a fixed-width fourth-order Runge-Kutta method.
#
# Args:
#   f: An array of functions defining the system of first-order ODEs such that
#       f[i](t, x_n) = x[i]_{n+1}
#   x0: The initial point.
#   t0: The initial time.
#   h: The step size.
#   n: The number of iterations.
#
# Throws:
#   ArgumentError: If the length of f is not equal to the length of x0.
#   ArgumentError: If n is not greater than zero.
#
# Returns: The vector of times, and the vector of points produced by the solver.
function rk4f{T<:AbstractFloat}(f::Array{Function,1}, x0::Array{T,1},
    t0::AbstractFloat, h::AbstractFloat, n::Integer)

    #Validate that the lengths of f and x0 are the same.
    if length(f) != length(x0)
        throw(ArgumentError("There must be one function per element of x0."));
    end

    #Validate that n is greater than zero.
    if n <= 0
        throw(ArgumentError("Number of iterations must be positive."));
    end

    #Validate that the step size is greater than zero.
    if h <= 0.0
        throw(ArgumentError("Step size must be positive."));
    end

    x = zeros(length(x0), n+1);
    t = zeros(1, n+1);

    #The first points are the initial conditions.
    x[:,1] = x0;
    t[1] = t0;

    for ii in 2:(n+1)
        t[ii], x[:,ii] = nextPoint(f, x[:,ii-1], t[ii-1], h);
    end

    return t,x
end #Close rk4f.

export rk4f

end
