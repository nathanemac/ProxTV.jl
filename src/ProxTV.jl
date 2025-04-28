module ProxTV

using OpenBLAS32_jll
using LinearAlgebra
using ShiftedProximalOperators

function __init__()
  config = LinearAlgebra.BLAS.lbt_get_config()
  if !any(lib -> lib.interface == :lp64, config.loaded_libs)
    LinearAlgebra.BLAS.lbt_forward(OpenBLAS32_jll.libopenblas_path)
  end
end

# import functions that we extend from ShiftedProximalOperators
import ShiftedProximalOperators.shift!
import ShiftedProximalOperators.shifted
import ShiftedProximalOperators.prox!

# export our new functions
export ProxTVContext
export InexactShiftedProximableFunction
export NormLp, ShiftedNormLp, NormTVp, ShiftedNormTVp
export prox!, shifted, shift!, TV, TVp_norm
export fun_name, fun_expr, fun_params

# main library functions
include("libproxtv.jl")
include("proxtv_utils.jl")

end
