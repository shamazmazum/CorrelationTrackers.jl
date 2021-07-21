"""
    Generic type for correlation function tracker.

See also: [`L2Tracker`](@ref), [`S2Tracker`](@ref),
[`SSTracker`](@ref).
"""
abstract type AbstractTracker{T} end

"""
    L2Tracker(phase)

Descriptor for line-segment correlation function for the phase
`phase`.

See also: [`S2Tracker`](@ref), [`SSTracker`](@ref),
[`AbstractTracker`](@ref).
"""
struct L2Tracker{T} <: AbstractTracker{T}
    phase :: T
end

"""
    S2Tracker(phase)

Descriptor for two-point correlation function for the phase
`phase`.

See also: [`L2Tracker`](@ref), [`SSTracker`](@ref),
[`AbstractTracker`](@ref).
"""
struct S2Tracker{T} <: AbstractTracker{T}
    phase :: T
end

"""
    SSTracker(phase)

Descriptor for surface-surface correlation function for the phase
`phase`.

See also: [`S2Tracker`](@ref), [`L2Tracker`](@ref),
[`AbstractTracker`](@ref).
"""
struct SSTracker{T} <: AbstractTracker{T}
    phase :: T
end

"""
    SVTracker(phase)

Descriptor for surface-void correlation function for the phase
`phase`.

See also: [`S2Tracker`](@ref), [`L2Tracker`](@ref),
[`SSTracker`](@ref), [`AbstractTracker`](@ref).
"""
struct SVTracker{T} <: AbstractTracker{T}
    phase :: T
end

# Utility functions
maybe_call_with_plans(slice :: AbstractArray{T},
                      data  :: S2Tracker{T};
                      plans :: Directional.S2FTPlans,
                      kwargs...) where T =
                          data(slice; plans = plans, kwargs...)
maybe_call_with_plans(slice :: AbstractArray{T},
                      data  :: AbstractTracker{T};
                      plans :: Directional.S2FTPlans,
                      kwargs...) where T =
                          data(slice; kwargs...)

# Is gradient update needed?
update_gradient_p(:: SSTracker)       = true
update_gradient_p(:: SVTracker)       = true
update_gradient_p(:: AbstractTracker) = false
