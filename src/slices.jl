# ? FIXME: Can this be merged with CorrelationFunctions.jl somehow?

# 2D
get_slice(a         :: AbstractArray{T, N},
          periodic  :: Bool,
          idx       :: NTuple{N, Int},
          direction :: Symbol) where {T, N} =
              get_slice(a, periodic, idx, Val(direction))

# Helper function for diagonals
function get_slice(a        :: AbstractArray{T, N},
                   periodic :: Bool,
                   idx      :: NTuple{N, Int},
                   forward  :: NTuple{N, Bool}) where {T, N}
    s = size(a)
    directions = (f ? 1 : -1 for f in forward)
    a = CircularArray(a)

    if periodic
        m = (length(a) == 3) ? idx[3] : idx[2]
        start = (f ? x - m + 1 : x + m - 1 for (f, x) in zip(forward, idx))
        indices = take(zip((countfrom(s, inc) for (s, inc) in zip(start, directions))...), s[1])
    else
        m = minimum(f ? x : s - x + 1 for (f, x, s) in zip(forward, idx, s))
        start = (f ? x - m + 1 : x + m - 1 for (f, x) in zip(forward, idx))
        indices = takewhile(x -> checkbounds(Bool, a, x...),
                            zip((countfrom(s, inc) for (s, inc) in zip(start, directions))...))
    end
    return [@inbounds a[index...] for index in indices], m
end

macro def_slicer(ndims, direction, expr)
    :(get_slice(a        :: AbstractArray{T, $ndims},
                periodic :: Bool,
                idx      :: NTuple{$ndims, Int},
                         :: Val{$direction}) where T = $expr)
end

# 2D
@def_slicer 2 :x let (x, y) = idx; a[:, y], x end
@def_slicer 2 :y let (x, y) = idx; a[x, :], y end
@def_slicer 2 :xy_main get_slice(a, periodic, idx, (true, true))
@def_slicer 2 :xy_anti get_slice(a, periodic, idx, (false, true))

# 3D
@def_slicer 3 :x let (x, y, z) = idx; a[:, y, z], x end
@def_slicer 3 :y let (x, y, z) = idx; a[x, :, z], y end
@def_slicer 3 :z let (x, y, z) = idx; a[x, y, :], z end

@def_slicer 3 :xy_main let (x, y, z) = idx; get_slice(a[:,:,z], periodic, (x, y), :xy_main) end
@def_slicer 3 :xy_anti let (x, y, z) = idx; get_slice(a[:,:,z], periodic, (x, y), :xy_anti) end
@def_slicer 3 :xz_main let (x, y, z) = idx; get_slice(a[:,y,:], periodic, (x, z), :xy_main) end
@def_slicer 3 :xz_anti let (x, y, z) = idx; get_slice(a[:,y,:], periodic, (x, z), :xy_anti) end
@def_slicer 3 :yz_main let (x, y, z) = idx; get_slice(a[x,:,:], periodic, (y, z), :xy_main) end
@def_slicer 3 :yz_anti let (x, y, z) = idx; get_slice(a[x,:,:], periodic, (y, z), :xy_anti) end

@def_slicer 3 :diag1 get_slice(a, periodic, idx, (true,  true,  true))
@def_slicer 3 :diag2 get_slice(a, periodic, idx, (false, true,  true))
@def_slicer 3 :diag3 get_slice(a, periodic, idx, (true,  false, true))
@def_slicer 3 :diag4 get_slice(a, periodic, idx, (true,  true,  false))
