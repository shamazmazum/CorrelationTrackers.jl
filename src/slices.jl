# ? FIXME: Can this be merged with CorrelationFunctions.jl somehow?

# 2D
get_slice(a         :: AbstractArray{T, N},
          idx       :: NTuple{N, Int},
          direction :: Symbol,) where {T, N} =
              get_slice(a, idx, Val(direction))

# Helper function for diagonals
function get_slice(a       :: AbstractArray{T, N},
                   idx     :: NTuple{N, Int},
                   forward :: NTuple{N, Bool}) where {T, N}
    s = size(a)
    m = minimum(f ? x : s - x + 1 for (f, x, s) in zip(forward, idx, s))
    directions = (f ? 1 : -1 for f in forward)
    start = (f ? x - m + 1 : x + m - 1 for (f, x) in zip(forward, idx))
    indices = takewhile(x -> checkbounds(Bool, a, x...),
        zip((countfrom(s, inc) for (s, inc) in zip(start, directions))...))
    return [@inbounds a[index...] for index in indices], m
end

get_slice(a :: AbstractArray{T, 2}, idx :: NTuple{2, Int}, :: Val{:x}) where T =
    let (x, y) = idx; a[:, y], x end

get_slice(a :: AbstractArray{T, 2}, idx :: NTuple{2, Int}, :: Val{:y}) where T =
    let (x, y) = idx; a[x, :], y end

get_slice(a :: AbstractArray{T, 2}, idx :: NTuple{2, Int}, :: Val{:xy_main}) where T =
    get_slice(a, idx, (true, true))

get_slice(a :: AbstractArray{T, 2}, idx :: NTuple{2, Int}, :: Val{:xy_anti}) where T =
    get_slice(a, idx, (false, true))

# 3D
get_slice(a :: AbstractArray{T, 3}, idx :: NTuple{3, Int}, :: Val{:x}) where T =
    let (x, y, z) = idx; a[:, y, z], x end

get_slice(a :: AbstractArray{T, 3}, idx :: NTuple{3, Int}, :: Val{:y}) where T =
    let (x, y, z) = idx; a[x, :, z], y end

get_slice(a :: AbstractArray{T, 3}, idx :: NTuple{3, Int}, :: Val{:z}) where T =
    let (x, y, z) = idx; a[x, y, :], z end

get_slice(a :: AbstractArray{T, 3}, idx :: NTuple{3, Int}, :: Val{:xy_main}) where T =
    let (x, y, z) = idx; get_slice(a[:,:,z], (x, y), :xy_main) end

get_slice(a :: AbstractArray{T, 3}, idx :: NTuple{3, Int}, :: Val{:xy_anti}) where T =
    let (x, y, z) = idx; get_slice(a[:,:,z], (x, y), :xy_anti) end

get_slice(a :: AbstractArray{T, 3}, idx :: NTuple{3, Int}, :: Val{:xz_main}) where T =
    let (x, y, z) = idx; get_slice(a[:,y,:], (x, z), :xy_main) end

get_slice(a :: AbstractArray{T, 3}, idx :: NTuple{3, Int}, :: Val{:xz_anti}) where T =
    let (x, y, z) = idx; get_slice(a[:,y,:], (x, z), :xy_anti) end

get_slice(a :: AbstractArray{T, 3}, idx :: NTuple{3, Int}, :: Val{:yz_main}) where T =
    let (x, y, z) = idx; get_slice(a[x,:,:], (y, z), :xy_main) end

get_slice(a :: AbstractArray{T, 3}, idx :: NTuple{3, Int}, :: Val{:yz_anti}) where T =
    let (x, y, z) = idx; get_slice(a[x,:,:], (y, z), :xy_anti) end

get_slice(a :: AbstractArray{T, 3}, idx :: NTuple{3, Int}, :: Val{:diag1}) where T =
    get_slice(a, idx, (true, true, true))

get_slice(a :: AbstractArray{T, 3}, idx :: NTuple{3, Int}, :: Val{:diag2}) where T =
    get_slice(a, idx, (false, true, true))

get_slice(a :: AbstractArray{T, 3}, idx :: NTuple{3, Int}, :: Val{:diag3}) where T =
    get_slice(a, idx, (true, false, true))

get_slice(a :: AbstractArray{T, 3}, idx :: NTuple{3, Int}, :: Val{:diag4}) where T =
    get_slice(a, idx, (true, true, false))
