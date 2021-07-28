var documenterSearchIndex = {"docs":
[{"location":"index.html#CorrelationTrackers.jl","page":"CorrelationTrackers.jl","title":"CorrelationTrackers.jl","text":"","category":"section"},{"location":"index.html","page":"CorrelationTrackers.jl","title":"CorrelationTrackers.jl","text":"CorrelationTrackers.jl package has means to do fast updates of correlation functions calculated by CorrelationFunctions.Directional module in CorrelationFunctions.jl. Correlation functions are recalculated when you change an element of the underlying array. Currently, only Directional.l2, Directional.s2 and Directional.surfsurf functions are supported.","category":"page"},{"location":"index.html","page":"CorrelationTrackers.jl","title":"CorrelationTrackers.jl","text":"An update after one element of the underlying system was changed is about 20000 times faster compared to a full recalculation of the correlation functions for that system.","category":"page"},{"location":"index.html#Examples","page":"CorrelationTrackers.jl","title":"Examples","text":"","category":"section"},{"location":"index.html","page":"CorrelationTrackers.jl","title":"CorrelationTrackers.jl","text":"This is a step-by-step guide on how to use CorrelationTrackers.jl:","category":"page"},{"location":"index.html","page":"CorrelationTrackers.jl","title":"CorrelationTrackers.jl","text":"For the beginning, create your system you want to calculate correlation functions for. Let it be two-dimensional two-phase system with size 20x30.\nPick combinations of correlation function + phase which you want to track, e.g S_2^1(r), L_2^1(r) and L_2^0(r). A combination of function and phase is represented by AbstractTracker structures, e.g. S2Tracker(1) will track S_2^(1)(r) function.\nCreate CorrelationTracker structure and access it as an ordinary array.\nObtain correlation functions at any time using functions from CorrelationFunctions.Directional. module.","category":"page"},{"location":"index.html","page":"CorrelationTrackers.jl","title":"CorrelationTrackers.jl","text":"using StatsBase\nusing PrettyTables\nusing Random\nusing CorrelationFunctions\nusing CorrelationTrackers\n\n# Create our system\nsystem = rand(MersenneTwister(348), 0:1, (20, 30))\n\n# Suppose we want to track $S_2^1(r)$ and $L_2^0(r)$.\ntracking = [S2Tracker(1), L2Tracker(1)]\n\n# Create the tracker. It may require some time if your system is big.\ntracker = CorrelationTracker(system; tracking = tracking)\n\n# Do some phase swapping in both the original system and the tracker\nsystem[12, 25]  = 1 - system[12, 25]\nsystem[11, 20]  = 1 - system[11, 20]\ntracker[12, 25] = 1 - tracker[12, 25]\ntracker[11, 20] = 1 - tracker[11, 20]\n\npretty_table(stdout, hcat(Directional.s2(system,  1) |> mean,\n                          Directional.l2(system,  1) |> mean,\n                          Directional.s2(tracker, 1) |> mean,\n                          Directional.l2(tracker, 1) |> mean);\n                     header = [\"S2 recalculated\", \"L2 recalculated\",\n                               \"S2 tracked\",      \"L2 tracked\"])","category":"page"},{"location":"index.html#API","page":"CorrelationTrackers.jl","title":"API","text":"","category":"section"},{"location":"index.html","page":"CorrelationTrackers.jl","title":"CorrelationTrackers.jl","text":"AbstractTracker\nS2Tracker\nL2Tracker\nSSTracker\nSVTracker\ndefault_trackers\nCorrelationTracker\ntracked_data\ntracked_length\ntracked_directions\nupdate_corrfns!\nrollback!","category":"page"},{"location":"index.html#CorrelationTrackers.AbstractTracker","page":"CorrelationTrackers.jl","title":"CorrelationTrackers.AbstractTracker","text":"Generic type for correlation function tracker.\n\nSee also: L2Tracker, S2Tracker, SSTracker.\n\n\n\n\n\n","category":"type"},{"location":"index.html#CorrelationTrackers.S2Tracker","page":"CorrelationTrackers.jl","title":"CorrelationTrackers.S2Tracker","text":"S2Tracker(phase)\n\nDescriptor for two-point correlation function for the phase phase.\n\nSee also: L2Tracker, SSTracker, AbstractTracker.\n\n\n\n\n\n","category":"type"},{"location":"index.html#CorrelationTrackers.L2Tracker","page":"CorrelationTrackers.jl","title":"CorrelationTrackers.L2Tracker","text":"L2Tracker(phase)\n\nDescriptor for line-segment correlation function for the phase phase.\n\nSee also: S2Tracker, SSTracker, AbstractTracker.\n\n\n\n\n\n","category":"type"},{"location":"index.html#CorrelationTrackers.SSTracker","page":"CorrelationTrackers.jl","title":"CorrelationTrackers.SSTracker","text":"SSTracker(phase)\n\nDescriptor for surface-surface correlation function for the phase phase.\n\nSee also: S2Tracker, L2Tracker, AbstractTracker.\n\n\n\n\n\n","category":"type"},{"location":"index.html#CorrelationTrackers.SVTracker","page":"CorrelationTrackers.jl","title":"CorrelationTrackers.SVTracker","text":"SVTracker(phase)\n\nDescriptor for surface-void correlation function for the phase phase.\n\nSee also: S2Tracker, L2Tracker, SSTracker, AbstractTracker.\n\n\n\n\n\n","category":"type"},{"location":"index.html#CorrelationTrackers.default_trackers","page":"CorrelationTrackers.jl","title":"CorrelationTrackers.default_trackers","text":"default_trackers(T)\n\nConstruct a vector of correlation functions which are tracked by default (that is S_2^1(x), L_2^1(x) and L_2^0(x)). T is the type of x.\n\n\n\n\n\n","category":"function"},{"location":"index.html#CorrelationTrackers.CorrelationTracker","page":"CorrelationTrackers.jl","title":"CorrelationTrackers.CorrelationTracker","text":"CorrelationTracker(system   :: AbstractArray{T, N}; \n                   tracking = default_trackers(T),\n                   periodic = false[, directions][, kwargs...])\n\nCreate correlation functions tracker.\n\nCreate correlation tracker for the array system. tracking is a vector of AbstractTracker structures which specify correlation functions you wish to track. periodic and direction have the same meaning as in the most functions in CorrelationFunctions.jl package. Additional arguments such as len may be passed in kwargs.\n\nReturned tracker supports interface of AbstractArray (e.g. you can perform element-wise read and write operations).\n\nExamples\n\njulia> let\n       system = rand(MersenneTwister(35), 0:1, (30, 10))\n       tracker = CorrelationTracker(system)\n       end\n30×10 CorrelationTracker{Int64, 2, Matrix{Int64}}:\n 0  1  0  1  1  0  0  1  1  0\n 1  1  1  0  0  0  0  0  1  1\n 0  0  0  0  0  0  1  1  0  1\n 1  1  1  0  1  1  1  0  1  0\n 0  1  0  0  1  0  0  1  1  1\n 0  0  0  0  0  0  1  0  1  1\n 0  0  1  0  1  1  0  1  0  1\n 1  0  0  1  0  0  1  0  1  0\n 0  1  1  0  0  1  1  1  1  1\n 0  0  1  1  1  1  0  0  0  0\n 0  0  1  1  0  0  1  1  1  0\n 0  1  0  0  0  1  0  0  1  0\n 1  0  0  1  0  0  1  1  0  1\n 0  1  0  1  0  0  1  1  1  0\n 1  1  0  1  1  1  0  1  0  1\n 1  1  1  0  0  0  0  1  0  1\n 1  0  0  1  0  0  1  1  1  0\n 0  0  0  1  0  0  0  1  1  0\n 1  0  1  0  1  0  0  0  1  0\n 1  0  0  1  0  0  0  0  0  1\n 1  1  1  0  1  0  1  0  1  1\n 0  1  0  1  1  0  0  1  0  1\n 0  0  0  1  0  0  1  1  1  1\n 0  0  1  1  1  1  0  1  1  0\n 1  0  1  1  0  0  0  0  0  1\n 1  1  0  1  0  1  1  0  1  0\n 0  1  1  0  0  1  1  0  1  0\n 0  1  0  0  1  0  0  1  0  0\n 1  1  0  0  1  1  0  0  0  1\n 0  0  1  1  0  1  1  1  1  0\n\n\n\n\n\n","category":"type"},{"location":"index.html#CorrelationTrackers.tracked_data","page":"CorrelationTrackers.jl","title":"CorrelationTrackers.tracked_data","text":"tracked_data(x :: CorrelationTracker)\n\nReturn an iterator over correlation function descriptors which are tracked by the tracker.\n\n\n\n\n\n","category":"function"},{"location":"index.html#CorrelationTrackers.tracked_length","page":"CorrelationTrackers.jl","title":"CorrelationTrackers.tracked_length","text":"tracked_length(x :: CorrelationTracker)\n\nReturn maximal tracked correlation length\n\nExamples\n\njulia> tracked_length(CorrelationTracker{Int,2}(rand(0:1, (50, 100))))\n25\n\n\n\n\n\n","category":"function"},{"location":"index.html#CorrelationTrackers.tracked_directions","page":"CorrelationTrackers.jl","title":"CorrelationTrackers.tracked_directions","text":"tracked_directions(x :: CorrelationTracker)\n\nReturn directions along which correlation functions are tracked.\n\nExamples\n\njulia> tracked_directions(CorrelationTracker{Int,2}(rand(0:1, (50, 100))))\n2-element Vector{Symbol}:\n :x\n :y\n\n\n\n\n\n","category":"function"},{"location":"index.html#CorrelationTrackers.update_corrfns!","page":"CorrelationTrackers.jl","title":"CorrelationTrackers.update_corrfns!","text":"update_corrfns!(tracker, value, index)\n\nThis function is equivalent to writing tracker[index] = value with exception that it also returns a rollback handle which can fastly bring the tracker to the previous state by calling rollback!.\n\nSee also: rollback!.\n\n\n\n\n\n","category":"function"},{"location":"index.html#CorrelationTrackers.rollback!","page":"CorrelationTrackers.jl","title":"CorrelationTrackers.rollback!","text":"rollback!(tracker, token)\n\nBring the system to the state before update_corrfns! was called. token must be an object returned by update_corrfns!.\n\nSee also: update_corrfns!.\n\n\n\n\n\n","category":"function"},{"location":"index.html#Caveats","page":"CorrelationTrackers.jl","title":"Caveats","text":"","category":"section"},{"location":"index.html","page":"CorrelationTrackers.jl","title":"CorrelationTrackers.jl","text":"Currently Directional.l2 and Directional.s2 return internal structures of supplied tracker. Do not modify them. This behavior can be changed to a safer one in the future.","category":"page"},{"location":"index.html","page":"CorrelationTrackers.jl","title":"CorrelationTrackers.jl","text":"Phase argument to SSTracker and SVTracker is currently ignored. It should be OK as long as tracked system is two-phase.","category":"page"}]
}
