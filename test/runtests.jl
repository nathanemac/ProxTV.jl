using ProxTV
using Test
using LinearAlgebra
using ShiftedProximalOperators

include("../src/libproxtv.jl")
include("../src/proxtv_utils.jl")

function simple_callback(s_ptr::Ptr{Cdouble}, s_length::Csize_t, delta_k::Cdouble, ctx_ptr::Ptr{Cvoid})
  return Cint(0)
end

@testset "ProxTV.jl" begin
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


  ctx = AlgorithmContextCallback(dualGap=dualGap)

  callback_pointer = @cfunction(simple_callback, Cint, (Ptr{Cdouble}, Csize_t, Cdouble, Ptr{Cvoid}))
  @test PN_LPp(y, lambda, x, info, n, p, ws, positive, ctx, callback_pointer) == 1 # 1 is the expected return value of the function. This means that the function has been executed successfully.

  @test TV(y, lambda, x, info, n, p, ws, ctx, callback_pointer) == 1
end
