# ? FIXME: Can this be merged with CorrelationFunctions.jl somehow?

get_slice(a :: AbstractArray{T, N}, direction :: Symbol, idx :: NTuple{N, Int}) where {T, N} =
    get_slice(a, Val(direction), idx)

get_slice(a :: AbstractArray{T, 2}, :: Val{:x}, idx :: NTuple{2, Int}) where T =
    let x = idx[1]; y = idx[2]; a[:, y], x end

get_slice(a :: AbstractArray{T, 2}, :: Val{:y}, idx :: NTuple{2, Int}) where T =
    let x = idx[1]; y = idx[2]; a[x, :], y end

get_slice(a :: AbstractArray{T, 3}, :: Val{:x}, idx :: NTuple{3, Int}) where T =
    let x = idx[1]; y = idx[2]; z = idx[3]; a[:, y, z], x end

get_slice(a :: AbstractArray{T, 3}, :: Val{:y}, idx :: NTuple{3, Int}) where T =
    let x = idx[1]; y = idx[2]; z = idx[3]; a[x, :, z], y end

get_slice(a :: AbstractArray{T, 3}, :: Val{:z}, idx :: NTuple{3, Int}) where T =
    let x = idx[1]; y = idx[2]; z = idx[3]; a[x, y, :], z end
