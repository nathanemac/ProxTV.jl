"""
    InexactShiftedProximableFunction

Abstract type for inexact proximal functions with shifts.
Used for functions that require iterative computation of their proximal operator.
"""
abstract type InexactShiftedProximableFunction end

export default_proxTV_callback
export default_proxTV_callback_v2
export default_proxTV_callback_v3

export update_prox_context!

# This is a helper struct to store the gradient and the proximal term in a compact way.
# It is used to avoid memory allocations when calling the proximal callback.
"""
    ModelFunction{V,P}

Helper structure to store the gradient and proximal term in a compact way.
Used to avoid memory allocations when calling the proximal callback.

# Fields
- `∇f::V`: gradient of the function
- `ψ::P`: proximal term
"""
mutable struct ModelFunction{V,P}
  ∇f::V  # gradient
  ψ::P   # proximal term
end

# define function that creates the structure ModelFunction
"""
    ModelFunction(∇f::V, ψ::Function) where {V<:AbstractVector}

Constructor for ModelFunction that creates a structure with a gradient vector and a proximal function.
"""
function ModelFunction(∇f::V, ψ::Function) where {V<:AbstractVector}
  return ModelFunction{V,Function}(∇f, ψ)
end

"""
    (m::ModelFunction)(d)

Evaluate the model function at point d by computing the sum of:
1. The inner product between the gradient and d
2. The proximal term evaluated at d
"""
function (m::ModelFunction)(d)
  return dot(m.∇f, d) + m.ψ(d)
end

### C++ callback functions

"""
    default_proxTV_callback(s_ptr::Ptr{Cdouble}, s_length::Csize_t, delta_k::Cdouble, ctx_ptr::Ptr{Cvoid})::Cint

Default callback function for ProxTV algorithm. Implements the stopping criterion based on the ratio between
the duality gap and the model decrease.

# Arguments
- `s_ptr`: Pointer to the current solution
- `s_length`: Length of the solution vector
- `delta_k`: Current duality gap
- `ctx_ptr`: Pointer to the ProxTVContext object

# Returns
- `Int32(1)` if the stopping criterion is satisfied
- `Int32(0)` otherwise

# Note
The stopping criterion is: δₖ ≤ ((1-κξ)/κξ) * ξₖ
where ξₖ is the model decrease at iteration k
"""
function default_proxTV_callback(
  s_ptr::Ptr{Cdouble},
  s_length::Csize_t,
  delta_k::Cdouble,
  ctx_ptr::Ptr{Cvoid},
)::Cint

  context = unsafe_pointer_to_objref(ctx_ptr)::ProxTVContext
  context.s_k .= unsafe_wrap(Vector{Float64}, s_ptr, s_length)

  # In-place operation to avoid memory allocations
  @. context.s_k_unshifted = context.s_k - context.shift

  # Computations without allocations
  ξk =
    context.hk - Float64(context.mk(context.s_k_unshifted)) +
    max(1, abs(context.hk)) * 10 * eps()
  condition = delta_k ≤ (1.0 - context.κξ) / context.κξ * ξk
  return condition ? Int32(1) : Int32(0)
end

"""
    default_proxTV_callback_v2(s_ptr::Ptr{Cdouble}, s_length::Csize_t, delta_k::Cdouble, ctx_ptr::Ptr{Cvoid})::Cint

Alternative callback function for ProxTV algorithm. Uses a fixed duality gap threshold from the context
and ensures the model decrease is non-negative.

# Arguments
- `s_ptr`: Pointer to the current solution
- `s_length`: Length of the solution vector
- `delta_k`: Current duality gap
- `ctx_ptr`: Pointer to the ProxTVContext object

# Returns
- `Int32(1)` if the stopping criterion is satisfied
- `Int32(0)` otherwise

# Note
The stopping criterion is: (δₖ ≤ dualGap) && (ξₖ ≥ 0)
where dualGap is fixed in the context
"""
function default_proxTV_callback_v2(
  s_ptr::Ptr{Cdouble},
  s_length::Csize_t,
  delta_k::Cdouble,
  ctx_ptr::Ptr{Cvoid},
)::Cint
  s_k = unsafe_wrap(Vector{Float64}, s_ptr, s_length; own = false)
  context = unsafe_pointer_to_objref(ctx_ptr)::ProxTVContext

  # In-place operation to avoid memory allocations
  @. context.s_k_unshifted = s_k - context.shift

  # Computations without allocations
  ξk = context.hk - context.mk(context.s_k_unshifted) + max(1, abs(context.hk)) * 10 * eps()

  condition = (delta_k ≤ context.dualGap) && (ξk ≥ 0)

  return condition ? Int32(1) : Int32(0)
end

"""
    default_proxTV_callback_v3(s_ptr::Ptr{Cdouble}, s_length::Csize_t, delta_k::Cdouble, ctx_ptr::Ptr{Cvoid})::Cint

Advanced callback function for ProxTV algorithm. Dynamically updates the duality gap threshold
based on the model decrease while ensuring it remains non-negative.

# Arguments
- `s_ptr`: Pointer to the current solution
- `s_length`: Length of the solution vector
- `delta_k`: Current duality gap
- `ctx_ptr`: Pointer to the ProxTVContext object

# Returns
- `Int32(1)` if the stopping criterion is satisfied
- `Int32(0)` otherwise

# Note
The stopping criterion is: (δₖ ≤ dualGap) && (ξₖ ≥ 0)
where dualGap is dynamically updated as: min(dualGap, ((1-κξ)/κξ) * ξₖ)
"""
function default_proxTV_callback_v3(
  s_ptr::Ptr{Cdouble},
  s_length::Csize_t,
  delta_k::Cdouble,
  ctx_ptr::Ptr{Cvoid},
)::Cint
  s_k = unsafe_wrap(Vector{Float64}, s_ptr, s_length; own = false)
  context = unsafe_pointer_to_objref(ctx_ptr)::ProxTVContext

  # In-place operation to avoid memory allocations
  @. context.s_k_unshifted = s_k - context.shift

  # Computations without allocations
  ξk = context.hk - context.mk(context.s_k_unshifted) + max(1, abs(context.hk)) * 10 * eps()

  aux = (1 - context.κξ) / context.κξ * ξk

  if aux < context.dualGap && aux ≥ 0
    context.dualGap = aux
  end

  condition = (delta_k ≤ context.dualGap) && (ξk ≥ 0)

  return condition ? Int32(1) : Int32(0)
end

## Structure for callback function in iR2/iR2N
"""
    ProxTVContext

Structure for managing the context of ProxTV operations, including memory allocations
and algorithm parameters.

# Fields
- `hk::Float64`: current step size
- `mk::ModelFunction`: model function (gradient + proximal term)
- `κξ::Float64`: control parameter between 1/2 and 1
- `shift::Vector{Float64}`: shift vector
- `s_k_unshifted::Vector{Float64}`: current unshifted solution
- `dualGap::Float64`: target duality gap
- `prox_stats::Any`: statistics about iterations
- `callback_pointer::Ptr{Cvoid}`: pointer to the C callback function
- `info::Vector{Float64}`: info array (3 elements)
- `temp_x::Vector{Float64}`: temporary vector for computations
- `y_shifted::Vector{Float64}`: for shifted versions
- `s::Vector{Float64}`: to store s = x - xk - sj
"""
mutable struct ProxTVContext
  hk::Float64
  mk::ModelFunction
  κξ::Float64
  shift::Vector{Float64}
  s_k_unshifted::Vector{Float64}
  dualGap::Float64
  prox_stats::Vector{Int64}  # for total number of iterations in ir2n, ir2 and prox
  callback_pointer::Ptr{Cvoid}  # pointer to the C callback function

  # Allocations for prox!
  info::Vector{Float64}  # info array (3 elements)
  temp_x::Vector{Float64}  # temporary vector for computations
  y_shifted::Vector{Float64}  # for shifted versions
  s::Vector{Float64}  # to store s = x - xk - sj
  s_k::Vector{Float64}  # to store s_k in the callback

  # Constructeur interne
  function ProxTVContext(
    hk::Float64,
    mk::ModelFunction,
    κξ::Float64,
    shift::Vector{Float64},
    s_k_unshifted::Vector{Float64},
    dualGap::Float64,
    prox_stats::Vector{Int64},
    callback_pointer::Ptr{Cvoid},
    info::Vector{Float64},
    temp_x::Vector{Float64},
    y_shifted::Vector{Float64},
    s::Vector{Float64},
    s_k::Vector{Float64},
  )
    ctx = new(
      hk,
      mk,
      κξ,
      shift,
      s_k_unshifted,
      dualGap,
      prox_stats,
      callback_pointer,
      info,
      temp_x,
      y_shifted,
      s,
      s_k,
    )
    return ctx
  end
end

# Constructeur externe
function ProxTVContext(
  n::Int;
  κξ = 0.75,
  dualGap = 0.0,
  callback::Function = default_proxTV_callback,
)
  # Vérification des paramètres
  n <= 0 && throw(ArgumentError("number of variables must be positive"))
  (κξ <= 1 / 2 || κξ >= 1) && throw(ArgumentError("κξ must be strictly between 1/2 and 1"))
  dualGap < 0 && throw(ArgumentError("dualGap must be nonnegative"))

  shift = zeros(n)
  s_k_unshifted = zeros(n)
  hk = 0.0
  mk = ModelFunction(zeros(n), x -> x)
  info = zeros(Float64, 3)
  temp_x = zeros(Float64, n)
  y_shifted = zeros(Float64, n)
  s = zeros(Float64, n)
  prox_stats = zeros(Int64, 3)
  s_k = zeros(Float64, n)

  # Convert the Julia callback function to a C function pointer
  callback_pointer =
    @cfunction(default_proxTV_callback, Cint, (Ptr{Cdouble}, Csize_t, Cdouble, Ptr{Cvoid}))

  return ProxTVContext(
    hk,
    mk,
    κξ,
    shift,
    s_k_unshifted,
    dualGap,
    prox_stats,
    callback_pointer,
    info,
    temp_x,
    y_shifted,
    s,
    s_k,
  )
end

### NormLp and ShiftedNormLp Implementation

"""
    NormLp{T1,T2}

Structure representing the Lp norm with parameter p and scaling factor λ.

# Fields
- `λ::T1`: scaling factor (scalar or array)
- `p::T2`: norm parameter (≥ 1)
- `context::ProxTVContext`: context for proximal computations

# Constructor
    NormLp(λ::T1, p::T2, context::ProxTVContext) where {T1,T2}

# Exceptions
- `ArgumentError` if λ < 0 (scalar) or if any element of λ is negative (array)
- `ArgumentError` if p < 1
"""
mutable struct NormLp{T1,T2}
  λ::T1
  p::T2
  context::ProxTVContext

  function NormLp(λ::T1, p::T2, context::ProxTVContext) where {T1,T2}
    if λ isa Real
      λ < 0 && error("λ must be nonnegative")
    elseif λ isa AbstractArray
      eltype(λ) <: Real || error("Elements of λ must be real")
      any(λ .< 0) && error("All elements of λ must be nonnegative")
    else
      error("λ must be a real scalar or array")
    end

    p >= 1 || error("p must be greater than or equal to one")
    new{T1,T2}(λ, p, context)
  end
end

"""
    prox!(y, h::NormLp, q, ν)

Evaluates the proximity operator of an Lp norm.

# Arguments
- `y`: Array in which to store the result
- `h`: NormLp object
- `q`: Vector to which the proximity operator is applied
- `ν`: Scaling factor

# Note
The duality gap at the solution is guaranteed to be less than `h.context.dualGap`
"""
function prox!(y::AbstractArray, h::NormLp, q::AbstractArray, ν::Real)
  n = length(y)

  # Create a new workspace for this call
  ws = newWorkspace(n)
  if ws === nothing
    throw(ErrorException("Failed to allocate workspace"))
  end

  try
    info = h.context.info

    # Adjust lambda to account for ν (multiply λ by ν)
    lambda_scaled = h.λ * ν

    positive = Int32(all(v -> v >= 0, y) ? 1 : 0)

    # Use the callback from context
    PN_LPp(
      q,
      lambda_scaled,
      y,
      info,
      n,
      h.p,
      ws,
      positive,
      h.context,
      h.context.callback_pointer,
    )

    # add the number of iterations in prox to the context object
    h.context.prox_stats[3] += Int64(info[1])

    return y
  finally
    # Always free the workspace
    freeWorkspace(ws)
  end
end

"""
    (h::NormLp)(x::AbstractArray)

Evaluate the Lp norm at point x, scaled by λ.
"""
function (h::NormLp)(x::AbstractArray)
  return h.λ * LPnorm(x, length(x), h.p)
end

"""
    ShiftedNormLp{R,T,V0,V1,V2}

Structure representing a shifted Lp norm.

# Fields
- `h::NormLp{R,T}`: underlying Lp norm
- `xk::V0`: first shift
- `sj::V1`: second shift
- `sol::V2`: temporary solution
- `shifted_twice::Bool`: indicates if a second shift has been applied
- `xsy::V2`: temporary vector for computations

# Note
The context is accessible via `h.context`
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

Creates a ShiftedNormLp object with initial shift xk.

# Arguments
- `h`: Lp norm to shift
- `xk`: shift vector

# Returns
A new ShiftedNormLp object
"""
shifted(h::NormLp{R,T}, xk::AbstractVector{R}) where {R<:Real,T<:Real} =
  ShiftedNormLp(h, xk, zero(xk), false)

"""
    shifted(ψ::ShiftedNormLp, sj::AbstractVector)

Creates a new ShiftedNormLp object by adding a second shift sj.

# Arguments
- `ψ`: already shifted Lp norm
- `sj`: second shift vector

# Returns
A new ShiftedNormLp object with both shifts
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

"""
    fun_name(ψ::ShiftedNormLp)
    fun_expr(ψ::ShiftedNormLp)
    fun_params(ψ::ShiftedNormLp)

Utility functions to get a string representation of the shifted Lp norm.
"""
fun_name(ψ::ShiftedNormLp) = "shifted Lp norm"
fun_expr(ψ::ShiftedNormLp) = "t ↦ λ * ‖xk + sj + t‖ₚ"
fun_params(ψ::ShiftedNormLp) = "xk = $(ψ.xk)\n" * " "^14 * "sj = $(ψ.sj)"

"""
    (ψ::ShiftedNormLp)(y::AbstractVector)

Evaluate the shifted Lp norm at point y.
"""
function (ψ::ShiftedNormLp)(y::AbstractVector)
  @. ψ.xsy = ψ.xk + ψ.sj + y
  return ψ.h(ψ.xsy)
end

"""
    prox!(y, ψ::ShiftedNormLp, q, ν)

Evaluates the proximity operator of a shifted Lp norm.

# Arguments
- `y`: Array in which to store the result
- `ψ`: ShiftedNormLp object
- `q`: Vector to which the proximity operator is applied
- `ν`: Scaling factor

# Note
Uses the context from `ψ.h` for computations
"""
function prox!(y::AbstractArray, ψ::ShiftedNormLp, q::AbstractArray, ν::Real)
  n = length(y)

  # Create a new workspace for this call
  ws = newWorkspace(n)
  if ws === nothing
    throw(ErrorException("Failed to allocate workspace"))
  end

  try
    # Use pre-existing allocations from h's context
    context = ψ.h.context
    info = context.info
    x = context.temp_x
    y_shifted = context.y_shifted
    s = context.s

    # Compute y_shifted = xk + sj + q
    @. y_shifted = ψ.xk + ψ.sj + q

    # Adjust lambda to account for ν (multiply λ by ν)
    lambda_scaled = ψ.h.λ * ν

    positive = Int32(all(v -> v >= 0, y_shifted) ? 1 : 0)
    if ψ.h.p == 1
      PN_LP1(y_shifted, lambda_scaled, x, info, length(y))
    elseif ψ.h.p == 2
      PN_LP2(y_shifted, lambda_scaled, x, info, length(y))
    elseif ψ.h.p == Inf
      PN_LPi(y_shifted, lambda_scaled, x, info, length(y), ws)
    else
      PN_LPp(
        y_shifted,
        lambda_scaled,
        x,
        info,
        length(y),
        ψ.h.p,
        ws,
        positive,
        context,
        context.callback_pointer,
      )
    end

    # Compute s = x - xk - sj
    @. s = x - ψ.xk - ψ.sj

    # Store the result in y
    y .= s

    # add the number of iterations in prox to the context object
    context.prox_stats[3] += Int64(info[1])

    return y
  finally
    # Always free the workspace
    freeWorkspace(ws)
  end
end

### NormTVp and ShiftedNormTVp Implementation

"""
    NormTVp{T1,T2}

Structure representing the Total Variation (TV) norm with parameter p and scaling factor λ.

# Fields
- `λ::T1`: scaling factor (scalar or array)
- `p::T2`: norm parameter (≥ 1)
- `context::ProxTVContext`: context for proximal computations

# Constructor
    NormTVp(λ::T1, p::T2, context::ProxTVContext) where {T1,T2}

# Exceptions
- `ArgumentError` if λ < 0 (scalar) or if any element of λ is negative (array)
- `ArgumentError` if p < 1
"""
mutable struct NormTVp{T1,T2}
  λ::T1
  p::T2
  context::ProxTVContext

  function NormTVp(λ::T1, p::T2, context::ProxTVContext) where {T1,T2}
    if λ isa Real
      λ < 0 && error("λ must be nonnegative")
    elseif λ isa AbstractArray
      eltype(λ) <: Real || error("Elements of λ must be real")
      any(λ .< 0) && error("All elements of λ must be nonnegative")
    else
      error("λ must be a real scalar or array")
    end

    p >= 1 || error("p must be greater than or equal to one")
    new{T1,T2}(λ, p, context)
  end
end

"""
    TVp_norm(x::AbstractArray, p::Real)

Computes the TVp norm of vector x with parameter p.

# Arguments
- `x`: input vector
- `p`: norm parameter

# Returns
The TVp norm value: (∑ᵢ |xᵢ₊₁ - xᵢ|ᵖ)^(1/p)
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
    (h::NormTVp)(x::AbstractArray)

Evaluate the TVp norm at point x, scaled by λ.
"""
function (h::NormTVp)(x::AbstractArray)
  return h.λ * TVp_norm(x, h.p)
end

"""
    ShiftedNormTVp{R,T,V0,V1,V2}

Structure representing a shifted TV norm.

# Fields
- `h::NormTVp{R,T}`: underlying TV norm
- `xk::V0`: first shift
- `sj::V1`: second shift
- `sol::V2`: temporary solution
- `shifted_twice::Bool`: indicates if a second shift has been applied
- `xsy::V2`: temporary vector for computations

# Note
The context is accessible via `h.context`
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

Creates a ShiftedNormTVp object with initial shift xk.

# Arguments
- `h`: TV norm to shift
- `xk`: shift vector

# Returns
A new ShiftedNormTVp object
"""
shifted(h::NormTVp{R,T}, xk::AbstractVector{R}) where {R<:Real,T<:Real} =
  ShiftedNormTVp(h, xk, zero(xk), false)

"""
    shifted(ψ::ShiftedNormTVp, sj::AbstractVector)

Creates a new ShiftedNormTVp object by adding a second shift sj.

# Arguments
- `ψ`: already shifted TV norm
- `sj`: second shift vector

# Returns
A new ShiftedNormTVp object with both shifts
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

"""
    fun_name(ψ::ShiftedNormTVp)
    fun_expr(ψ::ShiftedNormTVp)
    fun_params(ψ::ShiftedNormTVp)

Utility functions to get a string representation of the shifted TV norm.
"""
fun_name(ψ::ShiftedNormTVp) = "shifted TVp norm"
fun_expr(ψ::ShiftedNormTVp) = "t ↦ λ * TVp(xk + sj + t)"
fun_params(ψ::ShiftedNormTVp) = "xk = $(ψ.xk)\n" * " "^14 * "sj = $(ψ.sj)"

"""
    (ψ::ShiftedNormTVp)(y::AbstractVector)

Evaluate the shifted TV norm at point y.
"""
function (ψ::ShiftedNormTVp)(y::AbstractVector)
  @. ψ.xsy = ψ.xk + ψ.sj + y
  return ψ.h(ψ.xsy)
end

"""
    prox!(y, ψ::ShiftedNormTVp, q, ν)

Evaluates the proximity operator of a shifted TV norm.

# Arguments
- `y`: Array in which to store the result
- `ψ`: ShiftedNormTVp object
- `q`: Vector to which the proximity operator is applied
- `ν`: Scaling factor

# Note
Uses the context from `ψ.h` for computations
"""
function prox!(y::AbstractArray, ψ::ShiftedNormTVp, q::AbstractArray, ν::Real)
  n = length(y)

  # Create a new workspace for this call
  ws = newWorkspace(n)
  if ws === nothing
    throw(ErrorException("Failed to allocate workspace"))
  end

  try
    # Use pre-existing allocations from h's context
    context = ψ.h.context
    info = context.info
    x = context.temp_x
    y_shifted = context.y_shifted
    s = context.s

    # Compute y_shifted = xk + sj + q
    @. y_shifted = ψ.xk + ψ.sj + q

    # Adjust lambda to account for ν
    lambda_scaled = ψ.h.λ * ν

    # Call the TV function from ProxTV package
    TV(
      y_shifted,
      lambda_scaled,
      x,
      info,
      length(y),
      ψ.h.p,
      ws,
      context,
      context.callback_pointer,
    )

    # Compute s = x - xk - sj
    @. s = x - ψ.xk - ψ.sj

    # Store the result in y
    y .= s

    # add the number of iterations in prox to the context object
    context.prox_stats[3] += Int64(info[1])

    return y
  finally
    # Always free the workspace
    freeWorkspace(ws)
  end
end

"""
    shift!(ψ::Union{ShiftedNormLp,ShiftedNormTVp}, shift::Vector)

Updates the shift of a ShiftedNormLp or ShiftedNormTVp object.

# Arguments
- `ψ`: object to update
- `shift`: new shift vector

# Note
Updates either `xk` or `sj` depending on the value of `shifted_twice`
"""
function shift!(ψ::Union{ShiftedNormLp,ShiftedNormTVp}, shift::Vector{R}) where {R<:Real}
  if ψ.shifted_twice
    ψ.sj .= shift
  else
    ψ.xk .= shift
  end
  return ψ
end

### general utility functions

"""
    shifted(h::Union{NormLp, NormTVp}, xk::AbstractVector)

Creates a shifted version of a norm.

# Arguments
- `h`: norm to shift (either NormLp or NormTVp)
- `xk`: shift vector

# Returns
- A ShiftedNormLp if h is NormLp
- A ShiftedNormTVp if h is NormTVp

# Throws
- `ArgumentError` if h is neither NormLp nor NormTVp
"""
function shifted(h::Union{NormLp,NormTVp}, xk::AbstractVector)
  if h isa NormLp
    return ShiftedNormLp(h, xk, zero(xk), false)
  elseif h isa NormTVp
    return ShiftedNormTVp(h, xk, zero(xk), false)
  else
    throw(ArgumentError("The function h must be either NormLp or NormTVp"))
  end
end

"""
    prox!(y, ψ::Union{InexactShiftedProximableFunction, ShiftedProximableFunction}, q, ν)

Evaluates the proximity operator of a shifted regularizer.

# Arguments
- `y`: Array in which to store the result
- `ψ`: Either a ShiftedProximableFunction or an InexactShiftedProximableFunction
- `q`: Vector to which the proximity operator is applied
- `ν`: Scaling factor

# Returns
The solution is stored in the input vector `y`, which is also returned

# Throws
- `ErrorException` if ψ is neither ShiftedProximableFunction nor InexactShiftedProximableFunction
"""
function prox!(
  y,
  ψ::Union{InexactShiftedProximableFunction,ShiftedProximableFunction},
  q,
  ν,
)
  if ψ isa ShiftedProximableFunction
    # Call to exact prox!()
    return prox!(y, ψ, q, ν)
  elseif ψ isa InexactShiftedProximableFunction
    # Call to inexact prox!()
    return prox!(y, ψ, q, ν)
  else
    error("Type $(typeof(ψ)) not supported")
  end
end

### Update context before prox! call

"""
    update_prox_context!(solver, stats, ψ)

Updates the context of an InexactShiftedProximableFunction object before calling prox!.

# Arguments
- `solver`: solver object
- `ψ`: InexactShiftedProximableFunction object
"""
function update_prox_context!(solver, stats, ψ::InexactShiftedProximableFunction)
  update_prox_context!(solver, stats, ψ, Val(typeof(ψ)))
end

"""
    update_prox_context!(solver, stats, ψ, T::Val{<:ShiftedNormLp})

Updates the context of a ShiftedNormLp object before calling prox!.

# Arguments
- `solver`: solver object
- `stats`: stats object
- `ψ`: ShiftedNormLp object
- `T`: Type of the object
"""
function update_prox_context!(solver, stats, ψ, T::Val{<:ShiftedNormLp})
  ψ.h.context.hk = stats.solver_specific[:nonsmooth_obj]
  ψ.h.context.mk.∇f = solver.∇fk
  ψ.h.context.mk.ψ = d -> ψ(d)  # Use the evaluation function of ψ instead of the object itself
  @. ψ.h.context.shift = ψ.xk + ψ.sj
end

"""
    update_prox_context!(solver, ψ, T::Val{<:ShiftedNormTVp})

Updates the context of a ShiftedNormTVp object before calling prox!.

# Arguments
- `solver`: solver object
- `ψ`: ShiftedNormTVp object
- `T`: Type of the object
"""
function update_prox_context!(solver, stats, ψ, T::Val{<:ShiftedNormTVp})
  ψ.h.context.hk = stats.solver_specific[:nonsmooth_obj]
  ψ.h.context.mk.∇f = solver.∇fk
  ψ.h.context.mk.ψ = d -> ψ(d)  # Use the evaluation function of ψ instead of the object itself
  @. ψ.h.context.shift = ψ.xk + ψ.sj
end
