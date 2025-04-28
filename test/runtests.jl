using ProxTV
using Test
using LinearAlgebra
using ShiftedProximalOperators

include("../src/libproxtv.jl")
include("../src/proxtv_utils.jl")

function simple_callback(
  s_ptr::Ptr{Cdouble},
  s_length::Csize_t,
  delta_k::Cdouble,
  ctx_ptr::Ptr{Cvoid},
)
  return Cint(0)
end

# Basic internal tests
@testset "Internal functions" begin
  # test on a simple function :
  n = 4
  x = rand(n)
  p = 2.0
  @test isapprox(LPnorm(x, n, p), norm(x), atol = 1e-5)

  ## test on more advanced functions :

  # test PN_LPp
  n = 3
  y = [1.0, -2.0, 3.0]
  lambda = 0.1
  x = zeros(n)
  info = zeros(Float64, 3)
  ws = newWorkspace(n)
  positive = Int32(0)
  dualGap = 1e-4

  ctx = ProxTVContext(n, :lp, 1.7)

  callback_pointer =
    @cfunction(simple_callback, Cint, (Ptr{Cdouble}, Csize_t, Cdouble, Ptr{Cvoid}))
  @test PN_LPp(y, lambda, x, info, n, p, ws, positive, ctx, callback_pointer) == 1 # 1 is the expected return value of the function. This means that the function has been executed successfully.

  @test TV(y, lambda, x, info, n, p, ws, ctx, callback_pointer) == 1
end

@testset "ProxTVContext" begin
  n = 1000
  p = 1.7
  λ = 0.1
  @test_nowarn ctx_lp = ProxTVContext(n, :lp, p)
  @test_nowarn ctx_tvp = ProxTVContext(n, :tvp, p)

  ctx_lp = ProxTVContext(n, :lp, p)
  ctx_tvp = ProxTVContext(n, :tvp, p)

  @test_nowarn NormLp(λ, p, ctx_lp) # even though this is not the proper way to initialize the norm, see after.
  @test_nowarn NormTVp(λ, p, ctx_tvp)
end

@testset "Norms and ShiftedNorms " begin
  n = 1000
  x = rand(n)
  p = 1.7
  λ = 0.1
  q = rand(n)

  @test_nowarn h_lp = NormLp(λ, p, n) # the proper way to initialize the norm. The context is automatically created inside the NormLp constructor.
  @test_nowarn h_tvp = NormTVp(λ, p, n) # the proper way to initialize the norm. The context is automatically created inside the NormTVp constructor.

  h_lp = NormLp(λ, p, n)
  h_tvp = NormTVp(λ, p, n)

  y = similar(x)
  @test_nowarn h_lp(y)
  @test_nowarn h_tvp(y)

  @test_nowarn ψ_lp = shifted(h_lp, x)
  @test_nowarn ψ_tvp = shifted(h_tvp, x)

  ψ_lp = shifted(h_lp, x)
  ψ_tvp = shifted(h_tvp, x)

  @test_nowarn ψ_lp(y)
  @test_nowarn ψ_tvp(y)

  @test_nowarn shift!(ψ_lp, q)
  @test_nowarn shift!(ψ_tvp, q)
end

@testset "Proximal Operators" begin
  n = 1000
  x = rand(n)
  p = 1.7
  λ = 0.1

  h_lp = NormLp(λ, p, n)
  h_tvp = NormTVp(λ, p, n)

  y = similar(x)
  ν = 1e-5
  q = rand(n)

  @test_nowarn prox!(y, h_lp, q, ν)
  @test_nowarn prox!(y, h_tvp, q, ν)

  ψ_lp = shifted(h_lp, x)
  ψ_tvp = shifted(h_tvp, x)

  @test_nowarn prox!(y, ψ_lp, q, ν)
  @test_nowarn prox!(y, ψ_tvp, q, ν)
end

@testset "Callback" begin
  n = 1000
  s = rand(n)
  p = 1.7
  λ = 0.1
  delta_k = 0.1
  x = rand(n)
  s_length = Csize_t(length(s))
  delta_k = Cdouble(0.1)
  s_ptr = pointer(s)

  h_lp = NormLp(λ, p, n)
  ψ_lp = shifted(h_lp, x)
  ctx_ptr = Ptr{Cvoid}(pointer_from_objref(ψ_lp.h.context))
  @test_nowarn default_proxTV_callback_Lp(s_ptr, s_length, delta_k, ctx_ptr)

  h_tvp = NormTVp(λ, p, n)
  ψ_tvp = shifted(h_tvp, x)
  ctx_ptr = Ptr{Cvoid}(pointer_from_objref(ψ_tvp.h.context))
  @test_nowarn default_proxTV_callback_TVp(s_ptr, s_length, delta_k, ctx_ptr)
end
