module ProxTV

using OpenBLAS32_jll
using LAPACK_jll
using LinearAlgebra
# Include wrappers and the main library
include("libproxtv.jl")
# include("gen/wrapper.jl")

config = LinearAlgebra.BLAS.lbt_get_config()
if !any(lib -> lib.interface == :lp64, config.loaded_libs)
  LinearAlgebra.BLAS.lbt_forward(OpenBLAS32_jll.libopenblas_path)
end

end
