using ProxTV
using Test
using LinearAlgebra

@testset "ProxTV.jl" begin
  # test on a basic function :
  n = 4
  x = rand(n)
  p = 2.0
  @test isapprox(ProxTV.LPnorm(x, n, p), norm(x), atol = 1e-5)

  # test on a more advanced function :
  lambda = 0.18
  info = []
  y = rand(n)
  x = zeros(n)
  ws = ProxTV.newWorkspace(n)
  @test ProxTV.TV(y, lambda, x, info, n, p, ws) == 1 # 1 is the expected return value of the function

end
