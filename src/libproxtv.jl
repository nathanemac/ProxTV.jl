# This file contains the Julia bindings for the ProxTV library.

using proxTV_jll

function LPnorm(x, n, p)
  @ccall libproxtv.LPnorm(x::Ptr{Float64}, n::Int32, p::Float64)::Float64
end

function PN_LP1(y, lambda, x, info, n)
  @ccall libproxtv.PN_LP1(
    y::Ptr{Float64},
    lambda::Float64,
    x::Ptr{Float64},
    info::Ptr{Float64},
    n::Int32,
  )::Int32
end

function PN_LP2(y, lambda, x, info, n)
  @ccall libproxtv.PN_LP2(
    y::Ptr{Float64},
    lambda::Float64,
    x::Ptr{Float64},
    info::Ptr{Float64},
    n::Int32,
  )::Int32
end

struct Workspace
  n::Int32
  d::Ptr{Ptr{Float64}}
  maxd::Int32
  nd::Int32
  i::Ptr{Ptr{Int32}}
  maxi::Int32
  ni::Int32
  in::Ptr{Float64}
  out::Ptr{Float64}
  warm::Cshort
  warmDual::Ptr{Float64}
  warmLambda::Float64
end

function PN_LPinf(y, lambda, x, info, n, ws)
  @ccall libproxtv.PN_LPinf(
    y::Ptr{Float64},
    lambda::Float64,
    x::Ptr{Float64},
    info::Ptr{Float64},
    n::Int32,
    ws::Ptr{Workspace},
  )::Int32
end

# original PN_LPp function
function PN_LPp(y, lambda, x, info, n, p, ws, positive, ctx, callback)
  objGap = ctx.dualGap
  ctx_ptr = Ptr{Cvoid}(pointer_from_objref(ctx))
  @ccall libproxtv.PN_LPp(
    y::Ptr{Float64},
    lambda::Float64,
    x::Ptr{Float64},
    info::Ptr{Float64},
    n::Int32,
    p::Float64,
    ws::Ptr{Workspace},
    positive::Int32,
    objGap::Float64,
    ctx_ptr::Ptr{Cvoid},
    callback::Ptr{Cvoid},
  )::Int32
end

# overloaded PN_LPp function with less inputs
function PN_LPp(y, lambda, x, p, objGap)
  n = length(y)                  # works for nD signals
  info = []                      # stores the information about the execution of the function
  ws = ProxTV.newWorkspace(n)    # define a workspace for memory management
  positive = all(x -> x >= 0, y) # 0 if false, 1 if true
  @ccall libproxtv.PN_LPp(
    y::Ptr{Float64},
    lambda::Float64,
    x::Ptr{Float64},
    info::Ptr{Float64},
    n::Int32,
    p::Float64,
    ws::Ptr{Workspace},
    positive::Int32,
    objGap::Float64,
  )::Int32
end

function PN_LPp_v2(y, lambda, x, info, n, p, ws, positive)
  @ccall libproxtv.PN_LPp_v2(
    y::Ptr{Float64},
    lambda::Float64,
    x::Ptr{Float64},
    info::Ptr{Float64},
    n::Int32,
    p::Float64,
    ws::Ptr{Workspace},
    positive::Int32,
  )::Int32
end

function LP1_project(y, lambda, x, n, ws)
  @ccall libproxtv.LP1_project(
    y::Ptr{Float64},
    lambda::Float64,
    x::Ptr{Float64},
    n::Int32,
    ws::Ptr{Workspace},
  )::Int32
end

function LPp_project(y, lambda, x, info, n, p, ws)
  @ccall libproxtv.LPp_project(
    y::Ptr{Float64},
    lambda::Float64,
    x::Ptr{Float64},
    info::Ptr{Float64},
    n::Int32,
    p::Float64,
    ws::Ptr{Workspace},
  )::Int32
end

function solveLinearLP(z, n, p, lambda, s)
  @ccall libproxtv.solveLinearLP(
    z::Ptr{Float64},
    n::Int32,
    p::Float64,
    lambda::Float64,
    s::Ptr{Float64},
  )::Cvoid
end

# original TV function
function TV(y, lambda, x, info, n, p, ws, ctx, callback)
  objGap = ctx.dualGap
  ctx_ptr = Ptr{Cvoid}(pointer_from_objref(ctx))
  @ccall libproxtv.TV(
    y::Ptr{Float64},
    lambda::Float64,
    x::Ptr{Float64},
    info::Ptr{Float64},
    n::Int32,
    p::Float64,
    ws::Ptr{Workspace},
    objGap::Float64,
    ctx_ptr::Ptr{Cvoid},
    callback::Ptr{Cvoid},
  )::Int32
end

# overloaded TV function with less inputs
function TV(y, lambda, x, p)
  n = length(y)                  # works for nD signals
  info = []                      # stores the information about the execution of the function
  ws = ProxTV.newWorkspace(n)    # define a workspace for memory management
  @ccall libproxtv.TV(
    y::Ptr{Float64},
    lambda::Float64,
    x::Ptr{Float64},
    info::Ptr{Float64},
    n::Int32,
    p::Float64,
    ws::Ptr{Workspace},
  )::Int32
end

function PN_TV1(y, lambda, x, info, n, sigma, ws)
  @ccall libproxtv.PN_TV1(
    y::Ptr{Float64},
    lambda::Float64,
    x::Ptr{Float64},
    info::Ptr{Float64},
    n::Int32,
    sigma::Float64,
    ws::Ptr{Workspace},
  )::Int32
end

function linearizedTautString_TV1(y, lambda, x, n)
  @ccall libproxtv.linearizedTautString_TV1(
    y::Ptr{Float64},
    lambda::Float64,
    x::Ptr{Float64},
    n::Int32,
  )::Int32
end

function classicTautString_TV1(signal, n, lam, prox)
  @ccall libproxtv.classicTautString_TV1(
    signal::Ptr{Float64},
    n::Int32,
    lam::Float64,
    prox::Ptr{Float64},
  )::Int32
end

function hybridTautString_TV1(y, n, lambda, x)
  @ccall libproxtv.hybridTautString_TV1(
    y::Ptr{Float64},
    n::Int32,
    lambda::Float64,
    x::Ptr{Float64},
  )::Cvoid
end

function hybridTautString_TV1_custom(y, n, lambda, x, backtracksexp)
  @ccall libproxtv.hybridTautString_TV1_custom(
    y::Ptr{Float64},
    n::Int32,
    lambda::Float64,
    x::Ptr{Float64},
    backtracksexp::Float64,
  )::Cvoid
end

function classicTautString_TV1_offset(signal, n, lam, prox, offset)
  @ccall libproxtv.classicTautString_TV1_offset(
    signal::Ptr{Float64},
    n::Int32,
    lam::Float64,
    prox::Ptr{Float64},
    offset::Float64,
  )::Int32
end

function SolveTVConvexQuadratic_a1_nw(n, b, w, solution)
  @ccall libproxtv.SolveTVConvexQuadratic_a1_nw(
    n::Int32,
    b::Ptr{Float64},
    w::Float64,
    solution::Ptr{Float64},
  )::Cvoid
end

function PN_TV1_Weighted(Y, W, X, info, n, sigma, ws)
  @ccall libproxtv.PN_TV1_Weighted(
    Y::Ptr{Float64},
    W::Ptr{Float64},
    X::Ptr{Float64},
    info::Ptr{Float64},
    n::Int32,
    sigma::Float64,
    ws::Ptr{Workspace},
  )::Int32
end

function tautString_TV1_Weighted(y, lambda, x, n)
  @ccall libproxtv.tautString_TV1_Weighted(
    y::Ptr{Float64},
    lambda::Ptr{Float64},
    x::Ptr{Float64},
    n::Int32,
  )::Int32
end

function PN_TV1_Trend2_Weighted(Y, W, X, info, n, sigma, ws)
  @ccall libproxtv.PN_TV1_Trend2_Weighted(
    Y::Ptr{Float64},
    W::Ptr{Float64},
    X::Ptr{Float64},
    info::Ptr{Float64},
    n::Int32,
    sigma::Float64,
    ws::Ptr{Workspace},
  )::Int32
end

function SolveTVConvexQuadratic_a1(n, b, w, solution)
  @ccall libproxtv.SolveTVConvexQuadratic_a1(
    n::Int32,
    b::Ptr{Float64},
    w::Ptr{Float64},
    solution::Ptr{Float64},
  )::Cvoid
end

function more_TV2(y, lambda, x, info, n)
  @ccall libproxtv.more_TV2(
    y::Ptr{Float64},
    lambda::Float64,
    x::Ptr{Float64},
    info::Ptr{Float64},
    n::Int32,
  )::Int32
end

function morePG_TV2(y, lambda, x, info, n, ws)
  @ccall libproxtv.morePG_TV2(
    y::Ptr{Float64},
    lambda::Float64,
    x::Ptr{Float64},
    info::Ptr{Float64},
    n::Int32,
    ws::Ptr{Workspace},
  )::Int32
end

function PG_TV2(y, lambda, x, info, n)
  @ccall libproxtv.PG_TV2(
    y::Ptr{Float64},
    lambda::Float64,
    x::Ptr{Float64},
    info::Ptr{Float64},
    n::Int32,
  )::Int32
end

function DR2L1W_TV(M, N, unary, W1, W2, s, nThreads, maxit, info)
  @ccall libproxtv.DR2L1W_TV(
    M::Csize_t,
    N::Csize_t,
    unary::Ptr{Float64},
    W1::Ptr{Float64},
    W2::Ptr{Float64},
    s::Ptr{Float64},
    nThreads::Int32,
    maxit::Int32,
    info::Ptr{Float64},
  )::Int32
end

function GP_TVp(y, lambda, x, info, n, p, ws)
  @ccall libproxtv.GP_TVp(
    y::Ptr{Float64},
    lambda::Float64,
    x::Ptr{Float64},
    info::Ptr{Float64},
    n::Int32,
    p::Float64,
    ws::Ptr{Workspace},
  )::Int32
end

function OGP_TVp(y, lambda, x, info, n, p, ws)
  @ccall libproxtv.OGP_TVp(
    y::Ptr{Float64},
    lambda::Float64,
    x::Ptr{Float64},
    info::Ptr{Float64},
    n::Int32,
    p::Float64,
    ws::Ptr{Workspace},
  )::Int32
end

function FISTA_TVp(y, lambda, x, info, n, p, ws)
  @ccall libproxtv.FISTA_TVp(
    y::Ptr{Float64},
    lambda::Float64,
    x::Ptr{Float64},
    info::Ptr{Float64},
    n::Int32,
    p::Float64,
    ws::Ptr{Workspace},
  )::Int32
end

function FW_TVp(y, lambda, x, info, n, p, ws)
  @ccall libproxtv.FW_TVp(
    y::Ptr{Float64},
    lambda::Float64,
    x::Ptr{Float64},
    info::Ptr{Float64},
    n::Int32,
    p::Float64,
    ws::Ptr{Workspace},
  )::Int32
end

function GPFW_TVp(y, lambda, x, info, n, p, ws; objGap = 1e-5)
  @ccall libproxtv.GPFW_TVp(
    y::Ptr{Float64},
    lambda::Float64,
    x::Ptr{Float64},
    info::Ptr{Float64},
    n::Int32,
    p::Float64,
    ws::Ptr{Workspace},
    objGap::Float64,
  )::Int32
end

function PD2_TV(y, lambdas, norms, dims, x, info, ns, nds, npen, ncores, maxIters)
  @ccall libproxtv.PD2_TV(
    y::Ptr{Float64},
    lambdas::Ptr{Float64},
    norms::Ptr{Float64},
    dims::Ptr{Float64},
    x::Ptr{Float64},
    info::Ptr{Float64},
    ns::Ptr{Int32},
    nds::Int32,
    npen::Int32,
    ncores::Int32,
    maxIters::Int32,
  )::Int32
end

function proxDykstraTV2DWeighted(y, W, x, ns, nds, info, ncores, maxIters)
  @ccall libproxtv.proxDykstraTV2DWeighted(
    y::Ptr{Float64},
    W::Ptr{Float64},
    x::Ptr{Float64},
    ns::Ptr{Int32},
    nds::Int32,
    info::Ptr{Float64},
    ncores::Int32,
    maxIters::Int32,
  )::Int32
end

function DR2_TV(M, N, unary, W1, W2, norm1, norm2, s, nThreads, maxit, info)
  @ccall libproxtv.DR2_TV(
    M::Csize_t,
    N::Csize_t,
    unary::Ptr{Float64},
    W1::Float64,
    W2::Float64,
    norm1::Float64,
    norm2::Float64,
    s::Ptr{Float64},
    nThreads::Int32,
    maxit::Int32,
    info::Ptr{Float64},
  )::Int32
end

function CondatChambollePock2_TV(M, N, Y, lambda, X, alg, maxit, info)
  @ccall libproxtv.CondatChambollePock2_TV(
    M::Csize_t,
    N::Csize_t,
    Y::Ptr{Float64},
    lambda::Float64,
    X::Ptr{Float64},
    alg::Cshort,
    maxit::Int32,
    info::Ptr{Float64},
  )::Int32
end

function Yang2_TV(M, N, Y, lambda, X, maxit, info)
  @ccall libproxtv.Yang2_TV(
    M::Csize_t,
    N::Csize_t,
    Y::Ptr{Float64},
    lambda::Float64,
    X::Ptr{Float64},
    maxit::Int32,
    info::Ptr{Float64},
  )::Int32
end

function Kolmogorov2_TV(M, N, Y, lambda, X, maxit, info)
  @ccall libproxtv.Kolmogorov2_TV(
    M::Csize_t,
    N::Csize_t,
    Y::Ptr{Float64},
    lambda::Float64,
    X::Ptr{Float64},
    maxit::Int32,
    info::Ptr{Float64},
  )::Int32
end

function Yang3_TV(M, N, O, Y, lambda, X, maxit, info)
  @ccall libproxtv.Yang3_TV(
    M::Csize_t,
    N::Csize_t,
    O::Csize_t,
    Y::Ptr{Float64},
    lambda::Float64,
    X::Ptr{Float64},
    maxit::Int32,
    info::Ptr{Float64},
  )::Int32
end

function PD_TV(y, lambdas, norms, dims, x, info, ns, nds, npen, ncores, maxIters)
  @ccall libproxtv.PD_TV(
    y::Ptr{Float64},
    lambdas::Ptr{Float64},
    norms::Ptr{Float64},
    dims::Ptr{Float64},
    x::Ptr{Float64},
    info::Ptr{Float64},
    ns::Ptr{Int32},
    nds::Int32,
    npen::Int32,
    ncores::Int32,
    maxIters::Int32,
  )::Int32
end

function PDR_TV(y, lambdas, norms, dims, x, info, ns, nds, npen, ncores, maxIters)
  @ccall libproxtv.PDR_TV(
    y::Ptr{Float64},
    lambdas::Ptr{Float64},
    norms::Ptr{Float64},
    dims::Ptr{Float64},
    x::Ptr{Float64},
    info::Ptr{Float64},
    ns::Ptr{Int32},
    nds::Int32,
    npen::Int32,
    ncores::Int32,
    maxIters::Int32,
  )::Int32
end

function TV1D_denoise(input, output, width, lambda)
  @ccall libproxtv.TV1D_denoise(
    input::Ptr{Float64},
    output::Ptr{Float64},
    width::Int32,
    lambda::Float64,
  )::Cvoid
end

function TV1D_denoise_tautstring(input, output, width, lambda)
  @ccall libproxtv.TV1D_denoise_tautstring(
    input::Ptr{Float64},
    output::Ptr{Float64},
    width::Int32,
    lambda::Float64,
  )::Cvoid
end

function dp(n, y, lam, beta)
  @ccall libproxtv.dp(n::Int32, y::Ptr{Float64}, lam::Float64, beta::Ptr{Float64})::Cvoid
end

function radialProjection(x, n, norm, lambda)
  @ccall libproxtv.radialProjection(
    x::Ptr{Float64},
    n::Int32,
    norm::Float64,
    lambda::Float64,
  )::Cvoid
end

function newWorkspace(n)
  @ccall libproxtv.newWorkspace(n::Int32)::Ptr{Workspace}
end

function resetWorkspace(ws)
  @ccall libproxtv.resetWorkspace(ws::Ptr{Workspace})::Cvoid
end

function getDoubleWorkspace(ws)
  @ccall libproxtv.getDoubleWorkspace(ws::Ptr{Workspace})::Ptr{Float64}
end

function getIntWorkspace(ws)
  @ccall libproxtv.getIntWorkspace(ws::Ptr{Workspace})::Ptr{Int32}
end

function freeWorkspace(ws)
  @ccall libproxtv.freeWorkspace(ws::Ptr{Workspace})::Cvoid
end

function newWorkspaces(n, p)
  @ccall libproxtv.newWorkspaces(n::Int32, p::Int32)::Ptr{Ptr{Workspace}}
end

function freeWorkspaces(wa, p)
  @ccall libproxtv.freeWorkspaces(wa::Ptr{Ptr{Workspace}}, p::Int32)::Cvoid
end

function compareDoublesReversed(v1, v2)
  @ccall libproxtv.compareDoublesReversed(v1::Ptr{Cvoid}, v2::Ptr{Cvoid})::Int32
end

const STOP_PNLP = 0

const STOP_GAP_PNLP = 1.0e-5

const MAX_ITERS_PNLP = 1000

const SIGMA_PNLP = 0.05

const EPSILON_PNLP = 1.0e-15

const MIN_GRD_HESSIAN_PNLP = 1.0e-15

const LPPROJ_PSMALL = 1.002

const LPPROJ_PLARGE = 100

const MIN_STEP_PNLP = 1.0e-10

const STOP_PN = 1.0e-6

const SIGMA = 0.05

const MAX_ITERS_PN = 100

const STOP_MS = 1.0e-5

const STOP_MSSUB = 1.0e-6

const MAX_ITERS_MS = 100

const STOP_TVLP = 1.0e-5

const MAX_ITERS_TVLP = 10000

const MAX_ITERS_TVLPFW = 1000000

const MAX_ITERS_TVLPGPFW = 1000000

const MAX_NOIMP_TVLP = 10

const OBJGAP_LPPROX_TVLP = 1.0e-15

const LAMBDA_STEPS_TVLP = 1

const LAMBDA_REDUCTION_TVLP = 1000.0

const STOP_STEP_TVLP_FW = 1.0e-15

const FW_CYCLES_TVLP = 10

const MAX_NOIMP_TVLP_GPFW = 10FW_CYCLES_TVLP

const MIN_IMP_TVLP = 1.0e-10

const ENFORCED_CYCLES_TVLP_GPFW = 10

const STOP_PD = 1.0e-6

const MAX_ITERS_PD = 35

const MAX_ITERS_CONDAT = 2500

const STOP_CONDAT = 0

const MAX_ITERS_KOLMOGOROV = 2500

const STOP_KOLMOGOROV = 0

const MAX_ITERS_DR = 35

const MAX_ITERS_YANG = 35

const lapack_int = Int32

const DEBUG_N = 10

const N_INFO = 3

const INFO_ITERS = 0

const INFO_GAP = 1

const INFO_RC = 2

const EPSILON = 1.0e-10

const RC_OK = 0

const RC_ITERS = 1

const RC_STUCK = 2

const RC_ERROR = 3

const WS_MAX_MEMORIES = 100
