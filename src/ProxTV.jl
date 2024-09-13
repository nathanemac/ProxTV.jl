module ProxTV

using OpenBLAS32_jll
using LinearAlgebra
# Include wrappers and the main library
include("libproxtv.jl")
# include("gen/wrapper.jl")

function configure_blas()
  config = LinearAlgebra.BLAS.lbt_get_config()
  if !any(lib -> lib.interface == :lp64, config.loaded_libs)
    @info "No LP64 BLAS found, forwarding to OpenBLAS32_jll"
    LinearAlgebra.BLAS.lbt_forward(OpenBLAS32_jll.libopenblas_path)
  else
    @info "BLAS library is already loaded correctly"
  end
end

configure_blas()

end
