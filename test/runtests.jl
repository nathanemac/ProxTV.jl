using ProxTV
using Test
using LinearAlgebra

@testset "ProxTV.jl" begin
  # test on a simple function :
  n = 4
  x = rand(n)
  p = 2.0
  @test isapprox(ProxTV.LPnorm(x, n, p), norm(x), atol = 1e-5)

  # test on more advanced functions :

  n = rand(10:100)
  lambda = 0.15
  y = rand(n, n)
  x = zeros(n, n)
  ws = ProxTV.newWorkspace(n * n)
  objGap = 1e-5
  p = 1.32 # inexact prox computation : no closed-form for p = 1.32.

  @test ProxTV.PN_LPp(y, lambda, x, p, objGap) == 1 # 1 is the expected return value of the function. This means that the function has been executed successfully.

  @test ProxTV.TV(y, lambda, x, p) == 1
end
