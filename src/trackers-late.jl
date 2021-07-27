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
    tracker.corrdata[SSTracker{T}(phase)]

# Make AbstractTracker callable for convenience
(tracked :: L2Tracker{T})(tracker :: CorrelationTracker{T}) where T =
    tracker.corrdata[tracked]
(tracked :: L2Tracker{T})(system :: AbstractArray{T}; kwargs...) where T =
    Directional.l2(system, tracked.phase; kwargs...)

(tracked :: S2Tracker{T})(tracker :: CorrelationTracker{T}) where T =
    tracker.corrdata[tracked]
(tracked :: S2Tracker{T})(system :: AbstractArray{T}; kwargs...) where T =
    Directional.s2(system, tracked.phase; kwargs...)

(tracked :: SSTracker{T})(tracker :: CorrelationTracker{T}) where T =
    tracker.corrdata[tracked]
(tracked :: SSTracker{T})(system :: AbstractArray{T}; kwargs...) where T =
    Directional.surfsurf(system, tracked.phase; kwargs...)

(tracked :: SVTracker{T})(tracker :: CorrelationTracker{T}) where T =
    tracker.corrdata[tracked]
(tracked :: SVTracker{T})(system :: AbstractArray{T}; kwargs...) where T =
    Directional.surfvoid(system, tracked.phase; kwargs...)
