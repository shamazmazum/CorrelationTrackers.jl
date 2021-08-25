# ! FIXME: It should be safe to return internal structures as long as
# noone is going to modify them. I see no such scenario.

@doc raw"""
    Directional.l2(x :: CorrelationTracker, phase)

Return $L_2^{\text{phase}}$ function for an underlying system of the
tracker `x`.
"""
Directional.l2(tracker :: CorrelationTracker{T}, phase) where T =
    tracker.corrdata[L2Tracker{T}(phase)]

@doc raw"""
    Directional.s2(x :: CorrelationTracker, phase)

Return $S_2^{\text{phase}}$ function for an underlying system of the
tracker `x`.
"""
Directional.s2(tracker :: CorrelationTracker{T}, phase) where T =
    tracker.corrdata[S2Tracker{T}(phase)]

@doc raw"""
    Directional.surfsurf(x :: CorrelationTracker, phase)

Return $F_{ss}_2^{\text{phase}}$ function for an underlying system of the
tracker `x`.
"""
Directional.surfsurf(tracker :: CorrelationTracker{T}, phase) where T =
    tracker.corrdata[SSTracker{T}(phase)]

@doc raw"""
    Directional.surfvoid(x :: CorrelationTracker, phase)

Return $F_{sv}_2^{\text{phase}}$ function for an underlying system of the
tracker `x`.
"""
Directional.surfvoid(tracker :: CorrelationTracker{T}, phase) where T =
    tracker.corrdata[SVTracker{T}(phase)]

const TrackerAlike{T, N} = Union{CorrelationTracker{T, N}, ExtrapolatedData{T, N}}

"""
    CorrelationTracker(system, like)

Create a `CorrelationTracker` on top of array `system` with parameters
matching those in the tracker `like`.
"""
function CorrelationTracker(system :: AbstractArray{T, N},
                            like   :: TrackerAlike{T, N}) where {T, N}
    return CorrelationTracker(copy(system);
                              tracking   = like |> tracked_data |> collect,
                              periodic   = like.periodic,
                              directions = tracked_directions(like),
                              len        = tracked_length(like))
end

# Make AbstractTracker callable for convenience
(tracked :: L2Tracker{T})(tracker :: TrackerAlike{T}) where T =
    tracker.corrdata[tracked]
(tracked :: S2Tracker{T})(tracker :: TrackerAlike{T}) where T =
    tracker.corrdata[tracked]
(tracked :: SSTracker{T})(tracker :: TrackerAlike{T}) where T =
    tracker.corrdata[tracked]
(tracked :: SVTracker{T})(tracker :: TrackerAlike{T}) where T =
    tracker.corrdata[tracked]
(tracked :: CCTracker{T})(tracker :: TrackerAlike{T}) where T =
    tracker.corrdata[tracked]

"""
    cross_correlation(tracker :: CorrelationTracker, phase1, phase2)

Return precalculated values of cross correlation function for phases
`phase1` and `phase2`.
"""
cross_correlation(tracker :: CorrelationTracker{T}, phase1, phase2) where T =
    tracker.corrdata[CCTracker{T}(phase1, phase2)]
