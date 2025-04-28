"""
    InexactShiftedProximableFunction

Abstract type for inexact proximal functions with shifts.
Used for functions that require iterative computation of their proximal operator.
"""
abstract type InexactShiftedProximableFunction end

export default_proxTV_callback_TVp
export default_proxTV_callback_Lp
export default_proxTV_callback_v2
export default_proxTV_callback_v3

export update_prox_context!

"""
    default_proxTV_callback_TVp(s_ptr::Ptr{Cdouble}, s_length::Csize_t, delta_k::Cdouble, ctx_ptr::Ptr{Cvoid})::Cint

Default callback function used in ProxTV.jl for the TVp norm. Implements the stopping criterion based on the ratio between
the current duality gap δₖ and the model decrease ξₖ:
```
δₖ ≤ ((1-κξ)/κξ) * ξₖ
```
where κξ is a parameter between 1/2 and 1 that controls the precision of the proximal operator.

# Arguments
- `s_ptr`: Pointer to the current solution
- `s_length`: Length of the solution vector
- `delta_k`: Current duality gap
- `ctx_ptr`: Pointer to the ProxTVContext object

# Returns
- `Int32(1)` if the stopping criterion is satisfied
- `Int32(0)` otherwise

# Note
This function is optimized for ProxTV.jl. It is not intended to be used directly.
If the user wants to use a different function, it should be implemented with care.
"""
function default_proxTV_callback_TVp(
  s_ptr::Ptr{Cdouble},
  s_length::Csize_t,
  delta_k::Cdouble,
  ctx_ptr::Ptr{Cvoid},
)::Cint

  context = unsafe_pointer_to_objref(ctx_ptr)::ProxTVContext{typeof(TVp_norm)}
  @inbounds for i = 1:s_length
    context.s_k[i] = unsafe_load(s_ptr, i)
  end

  @. context.s_k_unshifted = context.s_k - context.shift
  n = Int(s_length)

  ψ_val::Float64 = TVp_norm(context.s_k, n, context.p)
  ϕ_val::Float64 = dot(context.∇fk, context.s_k_unshifted)
  mks = ϕ_val + ψ_val

  hk = context.hk::Float64
  κξ = context.κξ::Float64
  ξk::Float64 = hk - mks
  ratio::Float64 = (1.0 - κξ) / κξ
  condition::Bool = delta_k ≤ ratio * ξk

  return condition ? Int32(1) : Int32(0)
end


"""
    default_proxTV_callback_Lp(s_ptr::Ptr{Cdouble}, s_length::Csize_t, delta_k::Cdouble, ctx_ptr::Ptr{Cvoid})::Cint

Default callback function used in ProxTV.jl for the Lp norm. Implements the stopping criterion based on the ratio between
the current duality gap δₖ and the model decrease ξₖ:
```
δₖ ≤ ((1-κξ)/κξ) * ξₖ
```
where κξ is a parameter between 1/2 and 1 that controls the precision of the proximal operator.

# Arguments
- `s_ptr`: Pointer to the current solution
- `s_length`: Length of the solution vector
- `delta_k`: Current duality gap
- `ctx_ptr`: Pointer to the ProxTVContext object

# Returns
- `Int32(1)` if the stopping criterion is satisfied
- `Int32(0)` otherwise

# Note
This function is optimized for ProxTV.jl. It is not intended to be used directly.
If the user wants to use a different function, it should be implemented with care.
"""
function default_proxTV_callback_Lp(
  s_ptr::Ptr{Cdouble},
  s_length::Csize_t,
  delta_k::Cdouble,
  ctx_ptr::Ptr{Cvoid},
)::Cint

  context = unsafe_pointer_to_objref(ctx_ptr)::ProxTVContext{typeof(LPnorm)}
  @inbounds for i = 1:s_length
    context.s_k[i] = unsafe_load(s_ptr, i)
  end

  @. context.s_k_unshifted = context.s_k - context.shift
  n = Int(s_length)

  ψ_val::Float64 = LPnorm(context.s_k, n, context.p::Float64)
  ϕ_val::Float64 = dot(context.∇fk, context.s_k_unshifted)
  mks = ϕ_val + ψ_val

  hk = context.hk::Float64
  κξ = context.κξ::Float64
  ξk::Float64 = hk - mks
  ratio::Float64 = (1.0 - κξ) / κξ
  condition::Bool = delta_k ≤ ratio * ξk

  return condition ? Int32(1) : Int32(0)

end

function default_proxTV_callback_v2(
  s_ptr::Ptr{Cdouble},
  s_length::Csize_t,
  delta_k::Cdouble,
  ctx_ptr::Ptr{Cvoid},
)::Cint

  @error "default_proxTV_callback_v2 is deprecated. Use default_proxTV_callback_Lp or default_proxTV_callback_TVp instead."
  # s_k = unsafe_wrap(Vector{Float64}, s_ptr, s_length; own = false)
  # context = unsafe_pointer_to_objref(ctx_ptr)::ProxTVContext

  # # In-place operation to avoid memory allocations
  # @. context.s_k_unshifted = s_k - context.shift

  # # Computations without allocations
  # ξk = context.hk - context.mk(context.s_k_unshifted) + max(1, abs(context.hk)) * 10 * eps()

  # condition = (delta_k ≤ context.dualGap) && (ξk ≥ 0)

  # return condition ? Int32(1) : Int32(0)
end


function default_proxTV_callback_v3(
  s_ptr::Ptr{Cdouble},
  s_length::Csize_t,
  delta_k::Cdouble,
  ctx_ptr::Ptr{Cvoid},
)::Cint

  @error "default_proxTV_callback_v3 is deprecated. Use default_proxTV_callback_Lp or default_proxTV_callback_TVp instead."
  # s_k = unsafe_wrap(Vector{Float64}, s_ptr, s_length; own = false)
  # context = unsafe_pointer_to_objref(ctx_ptr)::ProxTVContext

  # # In-place operation to avoid memory allocations
  # @. context.s_k_unshifted = s_k - context.shift

  # # Computations without allocations
  # ξk = context.hk - context.mk(context.s_k_unshifted) + max(1, abs(context.hk)) * 10 * eps()

  # aux = (1 - context.κξ) / context.κξ * ξk

  # if aux < context.dualGap && aux ≥ 0
  #   context.dualGap = aux
  # end

  # condition = (delta_k ≤ context.dualGap) && (ξk ≥ 0)

  # return condition ? Int32(1) : Int32(0)
end

## Structure for callback function in iR2/iR2N
"""
    ProxTVContext

Structure for managing the context of ProxTV operations, including memory allocations
and algorithm parameters.

# Fields
- `hk::Float64`: current step size
- `h_symb::Symbol`: symbol of h to evaluate in the callback
- `∇fk::Vector{Float64}`: gradient of the function at the current point
- `h_fun::F`: function handle for h
- `p::Real`: parameter of the norm
- `κξ::Float64`: control parameter for the stopping criterion
- `shift::Vector{Float64}`: shift vector
- `s_k_unshifted::Vector{Float64}`: current unshifted solution
- `dualGap::Float64`: target duality gap
- `prox_stats::Any`: statistics about iterations
- `callback_pointer::Ptr{Cvoid}`: pointer to the C callback function
- `info::Vector{Float64}`: info array (3 elements)
- `temp_x::Vector{Float64}`: temporary vector for computations
- `y_shifted::Vector{Float64}`: for shifted versions
- `s::Vector{Float64}`: to store s = x - xk - sj
- `s_k::Vector{Float64}`: to store s_k in the callback
"""
mutable struct ProxTVContext{F}
  hk::Float64
  h_symb::Symbol
  ∇fk::Vector{Float64}
  h_fun::F
  p::Real
  κξ::Float64
  shift::Vector{Float64}
  s_k_unshifted::Vector{Float64}
  dualGap::Float64
  prox_stats::Vector{Int64}
  callback_pointer::Ptr{Cvoid}
  info::Vector{Float64}
  temp_x::Vector{Float64}
  y_shifted::Vector{Float64}
  s::Vector{Float64}
  s_k::Vector{Float64}

  function ProxTVContext{F}(
    hk::Float64,
    h_symb::Symbol,
    ∇fk::Vector{Float64},
    h_fun::F,
    p::Real,
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
  ) where {F}
    ctx = new{F}(
      hk,
      h_symb,
      ∇fk,
      h_fun,
      p,
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

function ProxTVContext(n::Int, h_symb::Symbol, p::Real; κξ = 0.75, dualGap = 0.0)
  n <= 0 && throw(ArgumentError("number of variables must be positive"))
  (κξ <= 1 / 2 || κξ >= 1) && throw(ArgumentError("κξ must be strictly between 1/2 and 1"))
  dualGap < 0 && throw(ArgumentError("dualGap must be nonnegative"))
  p >= 1 || throw(ArgumentError("p must be greater than or equal to one"))
  shift = zeros(n)
  s_k_unshifted = zeros(n)
  hk = 0.0
  ∇fk = similar(shift)
  if h_symb == :lp
    h_fun = LPnorm
  elseif h_symb == :tvp
    h_fun = TVp_norm
  end
  info = zeros(Float64, 3)
  temp_x = zeros(Float64, n)
  y_shifted = zeros(Float64, n)
  s = zeros(Float64, n)
  prox_stats = zeros(Int64, 3)
  s_k = zeros(Float64, n)

  if h_symb == :tvp
    callback_pointer = @cfunction(
      default_proxTV_callback_TVp,
      Cint,
      (Ptr{Cdouble}, Csize_t, Cdouble, Ptr{Cvoid})
    )
  elseif h_symb == :lp
    callback_pointer = @cfunction(
      default_proxTV_callback_Lp,
      Cint,
      (Ptr{Cdouble}, Csize_t, Cdouble, Ptr{Cvoid})
    )
  else
    error("h_symb must be :tvp or :lp")
  end

  return ProxTVContext{typeof(h_fun)}(
    hk,
    h_symb,
    ∇fk,
    h_fun,
    p,
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


"""
    NormLp{T1,T2,C}

Structure representing the Lp norm with parameter p and scaling factor λ.

# Fields
- `λ::T1`: scaling factor (scalar or array)
- `p::T2`: norm parameter (≥ 1)
- `context::C`: context for proximal computations

# Constructor
    NormLp(λ::T1, p::T2, context::C) where {T1,T2,C}

# Exceptions
- `ArgumentError` if λ < 0 (scalar) or if any element of λ is negative (array)
- `ArgumentError` if p < 1

# Note
This is not the most efficient way to define the Lp norm.
Instead, see `NormLp(λ::T1, p::T2, n::Int)` to avoid defining the context.
"""
mutable struct NormLp{T1,T2,C}
  λ::T1
  p::T2
  context::C

  function NormLp(λ::T1, p::T2, context::C) where {T1,T2,C}
    if λ isa Real
      λ < 0 && error("λ must be nonnegative")
    elseif λ isa AbstractArray
      eltype(λ) <: Real || error("Elements of λ must be real")
      any(λ .< 0) && error("All elements of λ must be nonnegative")
    else
      error("λ must be a real scalar or array")
    end

    if p ≠ context.p
      error("p in NormLp must be equal to the context p")
    end
    if context.h_symb ≠ :lp
      error("h_symb in NormLp must be :lp")
    end

    p >= 1 || error("p must be greater than or equal to one")
    new{T1,T2,C}(λ, p, context)
  end
end


"""
    NormLp(λ::T1, p::T2, n::Int) where {T1,T2}

Construct an Lp norm structure with parameter p and scaling factor λ.
The context is automatically defined.
"""
function NormLp(λ::T1, p::T2, n::Int) where {T1,T2}
  context = ProxTVContext(n, :lp, p)
  return NormLp(λ, p, context)
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
The quality of the proximal operator depends on κξ and the callback function, see `default_proxTV_callback_Lp`.
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
    freeWorkspace(ws)
  end
end

"""
    (h::NormLp)(x::AbstractArray)

Evaluate the Lp norm at point x, scaled by λ.
"""
function (h::NormLp)(x::AbstractVector{Float64})
  return h.λ * LPnorm(x, length(x), h.p)
end

"""
    ShiftedNormLp{R,T,C,V0,V1,V2}

Structure representing a shifted Lp norm.

# Fields
- `h::NormLp{R,T,C}`: underlying Lp norm
- `xk::V0`: first shift
- `sj::V1`: second shift
- `sol::V2`: temporary solution
- `shifted_twice::Bool`: indicates if a second shift has been applied
- `xsy::V2`: temporary vector for computations
"""
mutable struct ShiftedNormLp{
  R<:Real,
  T<:Real,
  C<:ProxTVContext,
  V0<:AbstractVector{R},
  V1<:AbstractVector{R},
  V2<:AbstractVector{R},
} <: InexactShiftedProximableFunction
  h::NormLp{R,T,C}
  xk::V0
  sj::V1
  sol::V2
  shifted_twice::Bool
  xsy::V2

  function ShiftedNormLp(
    h::NormLp{R,T,C},
    xk::AbstractVector{R},
    sj::AbstractVector{R},
    shifted_twice::Bool,
  ) where {R<:Real,T<:Real,C<:ProxTVContext}
    sol = similar(xk)
    xsy = similar(xk)
    new{R,T,C,typeof(xk),typeof(sj),typeof(sol)}(h, xk, sj, sol, shifted_twice, xsy)
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
shifted(h::NormLp{R,T,C}, xk::AbstractVector{R}) where {R<:Real,T<:Real,C<:ProxTVContext} =
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
  ψ::ShiftedNormLp{R,T,C,V0,V1,V2},
  sj::AbstractVector{R},
) where {
  R<:Real,
  T<:Real,
  V0<:AbstractVector{R},
  V1<:AbstractVector{R},
  V2<:AbstractVector{R},
  C<:ProxTVContext,
} = ShiftedNormLp(ψ.h, ψ.xk, sj, true)

"""
    fun_name(ψ::ShiftedNormLp)
    fun_expr(ψ::ShiftedNormLp)
    fun_params(ψ::ShiftedNormLp)

Utility functions to get a string representation of the shifted Lp norm.
"""
fun_name(ψ::ShiftedNormLp) = "shifted Lp norm"
fun_expr(ψ::ShiftedNormLp) = "t ↦ λ * ‖xk + sj + t‖ₚ"
fun_params(ψ::ShiftedNormLp) =
  "xk = $(ψ.xk)\n" *
  " "^14 *
  "sj = $(ψ.sj)\n" *
  " "^14 *
  "λ = $(ψ.h.λ)\n" *
  " "^14 *
  "p = $(ψ.h.p)"

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
Uses the context from `ψ.h` for computations.
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
    freeWorkspace(ws)
  end
end


"""
    NormTVp{T1,T2,C}

Structure representing the Total Variation (TV) norm with parameter p and scaling factor λ.

# Fields
- `λ::T1`: scaling factor (scalar or array)
- `p::T2`: norm parameter (≥ 1)
- `context::C`: context for proximal computations

# Constructor
    NormTVp(λ::T1, p::T2, context::C) where {T1,T2,C}

# Exceptions
- `ArgumentError` if λ < 0 (scalar) or if any element of λ is negative (array)
- `ArgumentError` if p < 1

# Note
This is not the most efficient way to define the TVp norm.
Instead, see `NormTVp(λ::T1, p::T2, n::Int)` to avoid defining the context.
"""
mutable struct NormTVp{T1,T2,C}
  λ::T1
  p::T2
  context::C

  function NormTVp(λ::T1, p::T2, context::C) where {T1,T2,C}
    if λ isa Real
      λ < 0 && error("λ must be nonnegative")
    elseif λ isa AbstractArray
      eltype(λ) <: Real || error("Elements of λ must be real")
      any(λ .< 0) && error("All elements of λ must be nonnegative")
    else
      error("λ must be a real scalar or array")
    end

    if p ≠ context.p
      error("p in NormTVp must be equal to the context p")
    end

    if context.h_symb ≠ :tvp
      error("h_symb in NormTVp must be :tvp")
    end

    p >= 1 || error("p must be greater than or equal to one")
    new{T1,T2,C}(λ, p, context)
  end
end

"""
    NormTVp(λ::T1, p::T2, n::Int) where {T1,T2}

Construct a TVp norm structure with parameter p and scaling factor λ.
The context is automatically defined.
"""
function NormTVp(λ::T1, p::T2, n::Int) where {T1,T2}
  context = ProxTVContext(n, :tvp, p)
  return NormTVp(λ, p, context)
end

"""
    TVp_norm(x::AbstractVector{Float64}, p::Real)

Computes the TVp norm of vector x with parameter p.

# Arguments
- `x`: input vector
- `p`: norm parameter

# Returns
The TVp norm value: (∑ᵢ |xᵢ₊₁ - xᵢ|ᵖ)^(1/p)
"""
@inline function TVp_norm(x::AbstractVector{Float64}, p::Float64)
  n = length(x)
  return TVp_norm(x, n, p)
end

@inline function TVp_norm(x::AbstractVector{Float64}, n::Int, p::Float64)
  s = 0.0
  @inbounds @simd for i = 1:(n-1)
    s += abs(x[i+1] - x[i])^p
  end
  return s^(1 / p)
end

"""
    (h::NormTVp)(x::AbstractArray)

Evaluate the TVp norm at point x, scaled by λ.
"""
function (h::NormTVp)(x::AbstractArray)
  return h.λ * TVp_norm(x, h.p)
end

"""
    ShiftedNormTVp{R,T,C,V0,V1,V2}

Structure representing a shifted TVp norm.

# Fields
- `h::NormTVp{R,T,C}`: underlying TVp norm
- `xk::V0`: first shift
- `sj::V1`: second shift
- `sol::V2`: temporary solution
- `shifted_twice::Bool`: indicates if a second shift has been applied
- `xsy::V2`: temporary vector for computations
"""
mutable struct ShiftedNormTVp{
  R<:Real,
  T<:Real,
  C<:ProxTVContext,
  V0<:AbstractVector{R},
  V1<:AbstractVector{R},
  V2<:AbstractVector{R},
} <: InexactShiftedProximableFunction
  h::NormTVp{R,T,C}
  xk::V0
  sj::V1
  sol::V2
  shifted_twice::Bool
  xsy::V2

  function ShiftedNormTVp(
    h::NormTVp{R,T,C},
    xk::AbstractVector{R},
    sj::AbstractVector{R},
    shifted_twice::Bool,
  ) where {R<:Real,T<:Real,C<:ProxTVContext}
    sol = similar(xk)
    xsy = similar(xk)
    new{R,T,C,typeof(xk),typeof(sj),typeof(sol)}(h, xk, sj, sol, shifted_twice, xsy)
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
shifted(h::NormTVp{R,T,C}, xk::AbstractVector{R}) where {R<:Real,T<:Real,C<:ProxTVContext} =
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
  ψ::ShiftedNormTVp{R,T,C,V0,V1,V2},
  sj::AbstractVector{R},
) where {
  R<:Real,
  T<:Real,
  V0<:AbstractVector{R},
  V1<:AbstractVector{R},
  V2<:AbstractVector{R},
  C<:ProxTVContext,
} = ShiftedNormTVp(ψ.h, ψ.xk, sj, true)

"""
    fun_name(ψ::ShiftedNormTVp)
    fun_expr(ψ::ShiftedNormTVp)
    fun_params(ψ::ShiftedNormTVp)

Utility functions to get a string representation of the shifted TV norm.
"""
fun_name(ψ::ShiftedNormTVp) = "shifted TVp norm"
fun_expr(ψ::ShiftedNormTVp) = "t ↦ λ * TVp(xk + sj + t)"
fun_params(ψ::ShiftedNormTVp) =
  "xk = $(ψ.xk)\n" *
  " "^14 *
  "sj = $(ψ.sj)\n" *
  " "^14 *
  "λ = $(ψ.h.λ)\n" *
  " "^14 *
  "p = $(ψ.h.p)"

"""
    (ψ::ShiftedNormTVp)(y::AbstractVector)

Evaluate the shifted TV norm at point y.
"""
function (ψ::ShiftedNormTVp)(y::AbstractVector)
  @. ψ.xsy = ψ.xk + ψ.sj + y
  return ψ.h(ψ.xsy)
end

"""
    prox!(y::AbstractArray, h::NormTVp, q::AbstractArray, ν::Real)

Evaluates the proximity operator of a TV norm.

# Arguments
- `y`: Array in which to store the result
- `h`: NormTVp object
- `q`: Vector to which the proximity operator is applied
- `ν`: Scaling factor

# Note
Uses the context from `h` for computations.
"""
function prox!(y::AbstractArray, h::NormTVp, q::AbstractArray, ν::Real)
  n = length(y)

  # Workspace temporaire
  ws = newWorkspace(n)
  ws === nothing && error("Failed to allocate workspace")

  try
    info = h.context.info

    # λ mis à l’échelle par ν
    λ_scaled = h.λ * ν

    # Appel à l’opérateur TVp de ProxTV
    TV(
      q,              # entrée
      λ_scaled,       # pénalité
      y,              # sortie (in-place)
      info,
      n,
      h.p,            # paramètre p de la TVp
      ws,
      h.context,
      h.context.callback_pointer,
    )

    # Statistiques : ajout des itérations effectuées
    h.context.prox_stats[3] += Int64(info[1])

    return y
  finally
    freeWorkspace(ws)
  end
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


"""
    update_prox_context!(solver, stats, ψ::ShiftedProximableFunction)

No update is needed for ShiftedProximableFunction objects therefore this function is a no-op.

# Arguments
- `solver`: solver object
- `stats`: stats object
- `ψ`: ShiftedProximableFunction object
"""
function update_prox_context!(solver, stats, ψ::ShiftedProximableFunction)
  return
end

"""
    update_prox_context!(solver, stats, ψ::ShiftedNormLp)

Updates the context of a ShiftedNormLp object before calling prox!.

# Arguments
- `solver`: solver object
- `stats`: stats object
- `ψ`: ShiftedNormLp object
"""
function update_prox_context!(solver, stats, ψ::ShiftedNormLp)
  ψ.h.context.hk = stats.solver_specific[:nonsmooth_obj]
  copy!(ψ.h.context.∇fk, solver.∇fk)
  @. ψ.h.context.shift = ψ.xk + ψ.sj
  return
end

"""
    update_prox_context!(solver, stats, ψ::ShiftedNormTVp)

Updates the context of a ShiftedNormTVp object before calling prox!.

# Arguments
- `solver`: solver object
- `stats`: stats object
- `ψ`: ShiftedNormTVp object
"""
function update_prox_context!(solver, stats, ψ::ShiftedNormTVp)
  ψ.h.context.hk = stats.solver_specific[:nonsmooth_obj]
  copy!(ψ.h.context.∇fk, solver.∇fk)
  @. ψ.h.context.shift = ψ.xk + ψ.sj
  return
end
