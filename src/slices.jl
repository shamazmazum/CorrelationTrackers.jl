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
    return [@inbounds a[index...] for index in indices]
end

macro def_slicer(ndims, direction, expr)
    :(get_slice(a        :: AbstractArray{T, $ndims},
                periodic :: Bool,
                idx      :: NTuple{$ndims, Int},
                _        :: Val{$direction}) where T = $expr)
end

# 2D
@def_slicer 2 :x let (_, y) = idx; a[:, y] end
@def_slicer 2 :y let (x, _) = idx; a[x, :] end
@def_slicer 2 :xy get_slice(a, periodic, idx, (true, true))
@def_slicer 2 :yx get_slice(a, periodic, idx, (false, true))

# 3D
@def_slicer 3 :x let (_, y, z) = idx; a[:, y, z] end
@def_slicer 3 :y let (x, _, z) = idx; a[x, :, z] end
@def_slicer 3 :z let (x, y, _) = idx; a[x, y, :] end

@def_slicer 3 :xy let (x, y, z) = idx; get_slice(a[:,:,z], periodic, (x, y), :xy) end
@def_slicer 3 :yx let (x, y, z) = idx; get_slice(a[:,:,z], periodic, (x, y), :yx) end
@def_slicer 3 :xz let (x, y, z) = idx; get_slice(a[:,y,:], periodic, (x, z), :xy) end
@def_slicer 3 :zx let (x, y, z) = idx; get_slice(a[:,y,:], periodic, (x, z), :yx) end
@def_slicer 3 :yz let (x, y, z) = idx; get_slice(a[x,:,:], periodic, (y, z), :xy) end
@def_slicer 3 :zy let (x, y, z) = idx; get_slice(a[x,:,:], periodic, (y, z), :yx) end

@def_slicer 3 :xyz get_slice(a, periodic, idx, (true,  true,  true))
@def_slicer 3 :yxz get_slice(a, periodic, idx, (false, true,  true))
@def_slicer 3 :xzy get_slice(a, periodic, idx, (true,  false, true))
@def_slicer 3 :zyx get_slice(a, periodic, idx, (true,  true,  false))

macro def_same_line(ndims, direction, expr)
    :(same_line_p(idx1 :: CartesianIndex{$ndims},
                  idx2 :: CartesianIndex{$ndims},
                  _    :: Val{$direction}) = $expr)
end

same_line_p(idx1      :: CartesianIndex{N},
            idx2      :: CartesianIndex{N},
            direction :: Symbol) where N =
                same_line_p(idx1, idx2, Val(direction))

@def_same_line 2 :x  idx1[2] == idx2[2]
@def_same_line 2 :y  idx1[1] == idx2[1]
@def_same_line 2 :xy (idx2[2] - idx1[2]) == (idx2[1] - idx1[1])
@def_same_line 2 :yx (idx2[2] - idx1[2]) == (idx1[1] - idx2[1])

@def_same_line 3 :x  idx1[2] == idx2[2] && idx1[3] == idx2[3]
@def_same_line 3 :y  idx1[1] == idx2[1] && idx1[3] == idx2[3]
@def_same_line 3 :z  idx1[1] == idx2[1] && idx1[2] == idx2[2]
@def_same_line 3 :xy (idx2[2] - idx1[2]) == (idx2[1] - idx1[1]) && idx1[3] == idx2[3]
@def_same_line 3 :yx (idx2[2] - idx1[2]) == (idx1[1] - idx2[1]) && idx1[3] == idx2[3]
@def_same_line 3 :xz (idx2[3] - idx1[3]) == (idx2[1] - idx1[1]) && idx1[2] == idx2[2]
@def_same_line 3 :zx (idx2[3] - idx1[3]) == (idx1[1] - idx2[1]) && idx1[2] == idx2[2]
@def_same_line 3 :yz (idx2[3] - idx1[3]) == (idx2[2] - idx1[2]) && idx1[1] == idx2[1]
@def_same_line 3 :zy (idx2[3] - idx1[3]) == (idx1[2] - idx2[2]) && idx1[1] == idx2[1]
@def_same_line 3 :xyz (idx2[3] - idx1[3]) == (idx2[2] - idx1[2]) == (idx2[1] - idx1[1])
@def_same_line 3 :yxz (idx2[3] - idx1[3]) == (idx2[2] - idx1[2]) == (idx1[1] - idx2[1])
@def_same_line 3 :xzy (idx2[3] - idx1[3]) == (idx1[2] - idx2[2]) == (idx2[1] - idx1[1])
@def_same_line 3 :zyx (idx1[3] - idx2[3]) == (idx2[2] - idx1[2]) == (idx2[1] - idx1[1])
