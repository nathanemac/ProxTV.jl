abstract type InexactShiftedProximableFunction end


## Structure for callback function in iR2N
mutable struct AlgorithmContextCallback
  hk::Float64
  mk::Function
  mk1::Function
  κξ::Float64
  shift::AbstractVector{Float64}
  s_k_unshifted::Vector{Float64}
  dualGap::Float64
  prox_stats # for total number of iterations in ir2n, ir2 and prox.
  flag_projLp::Int
  iters_prox_projLp::Int
end
function AlgorithmContextCallback(;hk=0.0, mk = x -> x, mk1 = x -> x, κξ = 0.0, shift = zeros(0), s_k_unshifted = zeros(0), dualGap = 0.0, prox_stats = [0.0, [], []], flag_projLp = 0, iters_prox_projLp = 1)
  AlgorithmContextCallback(hk, mk, mk1, κξ, shift, s_k_unshifted, dualGap, prox_stats, flag_projLp, iters_prox_projLp)
end


### NormLp and ShiftedNormLp Implementation

"""
    NormLp(λ::Real or AbstractArray, p::Real)

Represents the Lp norm with parameter `p` and scaling factor `λ`.
"""
struct NormLp{T1,T2}
    λ::T1
    p::T2

    function NormLp(λ::T1, p::T2) where {T1,T2}
        if λ isa Real
            λ < 0 && error("λ must be nonnegative")
        elseif λ isa AbstractArray
            eltype(λ) <: Real || error("Elements of λ must be real")
            any(λ .< 0) && error("All elements of λ must be nonnegative")
        else
            error("λ must be a real scalar or array")
        end

        p >= 1 || error("p must be greater than or equal to one")
        new{T1,T2}(λ, p)
    end
end

"""
    prox!(y, h::NormLp, q, ν, ctx_ptr, callback)

Evaluates inexactly the proximity operator of a Lp norm object.
The duality gap at the solution is guaranteed to be less than `dualGap`.

Inputs:
    - `y`: Array in which to store the result.
    - `h`: NormLp object.
    - `q`: Vector to which the proximity operator is applied.
    - `ν`: Scaling factor.
    - `ctx_ptr`: Pointer to the context object.
    - `callback`: Pointer to the callback function.
"""
function prox!(
        y::AbstractArray,
        h::NormLp,
        q::AbstractArray,
        ν::Real,
        context::AlgorithmContextCallback,
        callback::Ptr{Cvoid};
)

    n = length(y)
    ws = newWorkspace(n)

    # Allocate info array (based on C++ code)
    info = zeros(Float64, 3)

    # Adjust lambda to account for ν (multiply λ by ν)
    lambda_scaled = h.λ * ν

    positive = Int32(all(v -> v >= 0, y) ? 1 : 0)

    PN_LPp(q, lambda_scaled, y, info, n, h.p, ws, positive, context, callback)

    freeWorkspace(ws)

    # add the number of iterations in prox to the context object
    push!(context.prox_stats[3], info[1])

    return y
end

# Allows NormLp objects to be called as functions
function (h::NormLp)(x::AbstractArray)
    return h.λ * LPnorm(x, length(x), h.p)
end

"""
    ShiftedNormLp

A mutable struct representing a shifted NormLp function.
"""
mutable struct ShiftedNormLp{
    R<:Real,
    T<:Real,
    V0<:AbstractVector{R},
    V1<:AbstractVector{R},
    V2<:AbstractVector{R},
} <: InexactShiftedProximableFunction
    h::NormLp{R,T}
    xk::V0
    sj::V1
    sol::V2
    shifted_twice::Bool
    xsy::V2

    function ShiftedNormLp(
        h::NormLp{R,T},
        xk::AbstractVector{R},
        sj::AbstractVector{R},
        shifted_twice::Bool,
    ) where {R<:Real,T<:Real}
        sol = similar(xk)
        xsy = similar(xk)
        new{R,T,typeof(xk),typeof(sj),typeof(sol)}(h, xk, sj, sol, shifted_twice, xsy)
    end
end

"""
    shifted(h::NormLp, xk::AbstractVector)

Creates a ShiftedNormLp object with initial shift `xk`.
"""
shifted(h::NormLp{R,T}, xk::AbstractVector{R}) where {R<:Real,T<:Real} =
    ShiftedNormLp(h, xk, zero(xk), false)

"""
    shifted(ψ::ShiftedNormLp, sj::AbstractVector)

Creates a ShiftedNormLp object by adding a second shift `sj`.
"""
shifted(
    ψ::ShiftedNormLp{R,T,V0,V1,V2},
    sj::AbstractVector{R},
) where {
    R<:Real,
    T<:Real,
    V0<:AbstractVector{R},
    V1<:AbstractVector{R},
    V2<:AbstractVector{R},
} = ShiftedNormLp(ψ.h, ψ.xk, sj, true)

# Functions to get the name, expression, and parameters of the function
fun_name(ψ::ShiftedNormLp) = "shifted Lp norm"
fun_expr(ψ::ShiftedNormLp) = "t ↦ λ * ‖xk + sj + t‖ₚ"
fun_params(ψ::ShiftedNormLp) = "xk = $(ψ.xk)\n" * " "^14 * "sj = $(ψ.sj)"

"""
    shift!(ψ::ShiftedNormLp, shift::AbstractVector)

Updates the shift of a ShiftedNormLp object.
"""
function shift!(ψ::ShiftedNormLp, shift::AbstractVector{R}) where {R<:Real}
    if ψ.shifted_twice
        ψ.sj .= shift
    else
        ψ.xk .= shift
    end
    return ψ
end

# Allows ShiftedNormLp objects to be called as functions
function (ψ::ShiftedNormLp)(y::AbstractVector)
    @. ψ.xsy = ψ.xk + ψ.sj + y
    return ψ.h(ψ.xsy)
end

"""
    prox!(y, ψ::ShiftedNormLp, q, ν, ctx_ptr, callback)

Evaluates inexactly the proximity operator of a shifted Lp norm.
The duality gap at the solution is guaranteed to be less than `dualGap`.

Inputs:
    - `y`: Array in which to store the result.
    - `ψ`: ShiftedNormLp object.
    - `q`: Vector to which the proximity operator is applied.
    - `ν`: Scaling factor.
    - `ctx_ptr`: Pointer to the context object.
    - `callback`: Pointer to the callback function.
"""
function prox!(
    y::AbstractArray,
    ψ::ShiftedNormLp,
    q::AbstractArray,
    ν::Real,
    context::AlgorithmContextCallback,
    callback::Ptr{Cvoid};
)
    n = length(y)
    ws = newWorkspace(n) # to avoid unexplained memory leaks

    # Allocate info array (based on C++ code)
    info = zeros(Float64, 3)

    # Compute y_shifted = xk + sj + q
    y_shifted = ψ.xk .+ ψ.sj .+ q

    # Adjust lambda to account for ν (multiply λ by ν)
    lambda_scaled = ψ.h.λ * ν

    # Allocate the x vector to store the intermediate solution
    x = zeros(n)

    positive = Int32(all(v -> v >= 0, y_shifted) ? 1 : 0)
    if ψ.h.p == 1
        PN_LP1(y_shifted, lambda_scaled, x, info, n)
    elseif ψ.h.p == 2
        PN_LP2(y_shifted, lambda_scaled, x, info, n)
    elseif ψ.h.p == Inf
        PN_LPi(y_shifted, lambda_scaled, x, info, n, ws)
    else
        PN_LPp(y_shifted, lambda_scaled, x, info, n, ψ.h.p, ws, positive, context, callback)
    end
    # Compute s = x - xk - sj
    s = x .- ψ.xk .- ψ.sj

    # Store the result in y
    y .= s

    freeWorkspace(ws)

    # add the number of iterations in prox to the context object
    push!(context.prox_stats[3], info[1])

    return y
end



### NormTVp and ShiftedNormTVp Implementation

"""
    NormTVp(λ::Real or AbstractArray, p::Real)

Represents the Total Variation (TV) norm with parameter `p` and scaling factor `λ`.
"""
struct NormTVp{T1,T2}
    λ::T1
    p::T2

    function NormTVp(λ::T1, p::T2) where {T1,T2}
        if λ isa Real
            λ < 0 && error("λ must be nonnegative")
        elseif λ isa AbstractArray
            eltype(λ) <: Real || error("Elements of λ must be real")
            any(λ .< 0) && error("All elements of λ must be nonnegative")
        else
            error("λ must be a real scalar or array")
        end

        p >= 1 || error("p must be greater than or equal to one")
        new{T1,T2}(λ, p)
    end
end

"""
    TVp_norm(x::AbstractArray, p::Real)

Computes the TVp norm of vector `x` with parameter `p`.
"""
function TVp_norm(x::AbstractArray, p::Real)
    n = length(x)
    tvp_sum = 0.0
    for i = 1:(n-1)
        tvp_sum += abs(x[i+1] - x[i])^p
    end
    return tvp_sum^(1 / p)
end

"""
    prox!(y, h::NormTVp, q, ν; dualGap=1e-5)

Evaluates inexactly the proximity operator of a TVp norm object.
The duality gap at the solution is guaranteed to be less than `dualGap`.

Inputs:
    - `y`: Array in which to store the result.
    - `h`: NormTVp object.
    - `q`: Vector to which the proximity operator is applied.
    - `ν`: Scaling factor.
    - `ctx_ptr`: Pointer to the context object.
    - `callback`: Pointer to the callback function.
"""
function prox!(
        y::AbstractArray,
        h::NormTVp,
        q::AbstractArray,
        ν::Real,
        context::AlgorithmContextCallback,
        callback::Ptr{Cvoid})

    n = length(y)
    ws = newWorkspace(n)

    # Allocate info array (based on C++ code)
    info = zeros(Float64, 3)

    # Adjust λ by ν
    lambda_scaled = h.λ * ν

    TV(q, lambda_scaled, y, info, n, h.p, ws, context, callback)

    freeWorkspace(ws)

    # add the number of iterations in prox to the context object
    push!(context.prox_stats[3], info[1])

    return y
end

# Allows NormTVp objects to be called as functions
function (h::NormTVp)(x::AbstractArray)
    return h.λ * TVp_norm(x, h.p)
end

"""
    ShiftedNormTVp

A mutable struct representing a shifted NormTVp function.
"""
mutable struct ShiftedNormTVp{
    R<:Real,
    T<:Real,
    V0<:AbstractVector{R},
    V1<:AbstractVector{R},
    V2<:AbstractVector{R},
} <: InexactShiftedProximableFunction
    h::NormTVp{R,T}
    xk::V0
    sj::V1
    sol::V2
    shifted_twice::Bool
    xsy::V2

    function ShiftedNormTVp(
        h::NormTVp{R,T},
        xk::AbstractVector{R},
        sj::AbstractVector{R},
        shifted_twice::Bool,
    ) where {R<:Real,T<:Real}
        sol = similar(xk)
        xsy = similar(xk)
        new{R,T,typeof(xk),typeof(sj),typeof(sol)}(h, xk, sj, sol, shifted_twice, xsy)
    end
end

"""
    shifted(h::NormTVp, xk::AbstractVector)

Creates a ShiftedNormTVp object with initial shift `xk`.
"""
shifted(h::NormTVp{R,T}, xk::AbstractVector{R}) where {R<:Real,T<:Real} =
    ShiftedNormTVp(h, xk, zero(xk), false)

"""
    shifted(ψ::ShiftedNormTVp, sj::AbstractVector)

Creates a ShiftedNormTVp object by adding a second shift `sj`.
"""
shifted(
    ψ::ShiftedNormTVp{R,T,V0,V1,V2},
    sj::AbstractVector{R},
) where {
    R<:Real,
    T<:Real,
    V0<:AbstractVector{R},
    V1<:AbstractVector{R},
    V2<:AbstractVector{R},
} = ShiftedNormTVp(ψ.h, ψ.xk, sj, true)

# Functions to get the name, expression, and parameters of the function
fun_name(ψ::ShiftedNormTVp) = "shifted TVp norm"
fun_expr(ψ::ShiftedNormTVp) = "t ↦ λ * TVp(xk + sj + t)"
fun_params(ψ::ShiftedNormTVp) = "xk = $(ψ.xk)\n" * " "^14 * "sj = $(ψ.sj)"

"""
    shift!(ψ::ShiftedNormTVp, shift::AbstractVector)

Updates the shift of a ShiftedNormTVp object.
"""
function shift!(ψ::ShiftedNormTVp, shift::AbstractVector{R}) where {R<:Real}
    if ψ.shifted_twice
        ψ.sj .= shift
    else
        ψ.xk .= shift
    end
    return ψ
end

# Allows ShiftedNormTVp objects to be called as functions
function (ψ::ShiftedNormTVp)(y::AbstractVector)
    @. ψ.xsy = ψ.xk + ψ.sj + y
    return ψ.h(ψ.xsy)
end

"""
    prox!(y, ψ::ShiftedNormTVp, q, σ, ctx_ptr, callback)

Evaluates inexactly the proximity operator of a shifted TVp norm.
The duality gap at the solution is guaranteed to be less than `dualGap`.

Inputs:
    - `y`: Array in which to store the result.
    - `ψ`: ShiftedNormTVp object.
    - `q`: Vector to which the proximity operator is applied.
    - `ν`: Scaling factor.
    - `ctx_ptr`: Pointer to the context object.
    - `callback`: Pointer to the callback function.
"""
function prox!(y::AbstractArray, ψ::ShiftedNormTVp, q::AbstractArray, ν::Real, context::AlgorithmContextCallback, callback::Ptr{Cvoid})
    n = length(y)
    ws = newWorkspace(n)

    # Allocate info array (based on C++ code)
    info = zeros(Float64, 3)

    # Compute y_shifted = xk + sj + q
    y_shifted = ψ.xk .+ ψ.sj .+ q

    # Adjust lambda to account for ν (multiply λ by ν)
    lambda_scaled = ψ.h.λ * ν

    # Allocate the x vector to store the intermediate solution
    x = similar(y)

    # Call the TV function from ProxTV package
    TV(y_shifted, lambda_scaled, x, info, n, ψ.h.p, ws, context, callback)

    # Compute s = x - xk - sj
    s = x .- ψ.xk .- ψ.sj

    # Store the result in y
    y .= s

    freeWorkspace(ws)

    # add the number of iterations in prox to the context object
    push!(context.prox_stats[3], info[1])

    return y
end



### general utility functions

"""
    shifted(h::Union{NormLp, NormTVp}, xk::AbstractVector)

Creates a shifted version of `h` depending on its type.
If `h` is of type `NormLp`, returns a `ShiftedNormLp`.
If `h` is of type `NormTVp`, returns a `ShiftedNormTVp`.
"""
function shifted(h::Union{NormLp, NormTVp}, xk::AbstractVector)
    if h isa NormLp
        return ShiftedNormLp(h, xk, zero(xk), false)
    elseif h isa NormTVp
        return ShiftedNormTVp(h, xk, zero(xk), false)
    else
        throw(ArgumentError("The function h must be either NormLp or NormTVp"))
    end
end

"""
    prox!(y, ψ::Union{InexactShiftedProximableFunction, ShiftedProximableFunction}, q, ν; ctx_ptr, callback)

Evaluates the proximity operator of a shifted regularizer, choosing between exact and inexact calculations based on the type of `ψ` and the presence of pointers.

- If `ψ` is a `ShiftedProximableFunction` and pointers are not provided, computes the **exact** proximity operator.
- If `ψ` is an `InexactShiftedProximableFunction` and pointers are provided, computes the **inexact** proximity operator by calling the callback.

Inputs:
    - `y`: Array in which to store the result.
    - `ψ`: Either a `ShiftedProximableFunction` (for exact prox) or an `InexactShiftedProximableFunction` (for inexact prox).
    - `q`: Vector to which the proximity operator is applied.
    - `ν`: Scaling factor.
    - `ctx_ptr`: Pointer to the context object.
    - `callback`: Pointer to the callback function.
Outputs:
    - The solution is stored in the input vector `y`, which is also returned.

Errors:
    - Raises an error if `ψ` is of type `ShiftedProximableFunction` and pointers are provided, or if `ψ` is of type `InexactShiftedProximableFunction` and pointers are not provided.
"""
function prox!(y, ψ::Union{InexactShiftedProximableFunction, ShiftedProximableFunction
  }, q, ν, ctx_ptr, callback)
    if ψ isa ShiftedProximableFunction
        # Call to exact prox!()
        return prox!(y, ψ, q, ν)
    elseif ψ isa InexactShiftedProximableFunction
        # Call to inexact prox!()
        return prox!(y, ψ, q, ν, ctx_ptr, callback)

    else
        error("Combination of ψ::$(typeof(ψ)) presence/lack of pointers is not a valid call to prox!.
        Please provide pointers for InexactShiftedProximableFunction or omit them for ShiftedProximableFunction.")
    end
end
