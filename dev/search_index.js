var documenterSearchIndex = {"docs":
[{"location":"index.html#CorrelationTrackers.jl","page":"CorrelationTrackers.jl","title":"CorrelationTrackers.jl","text":"","category":"section"},{"location":"index.html","page":"CorrelationTrackers.jl","title":"CorrelationTrackers.jl","text":"CorrelationTrackers.jl package has a means to do fast updates of correlation functions calculated by CorrelationFunctions.Directional module in CorrelationFunctions.jl. Correlation functions are recalculated when you change an element of the underlying array. Currently, only Directional.l2 and Directional.s2 functions are supported.","category":"page"},{"location":"index.html","page":"CorrelationTrackers.jl","title":"CorrelationTrackers.jl","text":"An update after one element of the underlying system was changed is about 20000 times faster compared to a full recalculation of the correlation functions for that system.","category":"page"},{"location":"index.html#Examples","page":"CorrelationTrackers.jl","title":"Examples","text":"","category":"section"},{"location":"index.html","page":"CorrelationTrackers.jl","title":"CorrelationTrackers.jl","text":"This is a step-by-step guide on how to use CorrelationTrackers.jl:","category":"page"},{"location":"index.html","page":"CorrelationTrackers.jl","title":"CorrelationTrackers.jl","text":"For the beginning, create your system you want to calculate correlation functions for. Let it be two-dimensional two-phase system with size 20x30.\nPick combinations of correlation function + phase which you want to track, e.g S_2^1(r), L_2^1(r) and L_2^0(r). A combination of function and phase is represented in TrackerData structure, e.g. TrackerData(Directional.s2, 1) will track S_2^1(r) function.\nCreate CorrelationTracker structure and access it as an ordinary array.\nObtain correlation functions at any time using functions from CorrelationFunctions.Directional. module.","category":"page"},{"location":"index.html","page":"CorrelationTrackers.jl","title":"CorrelationTrackers.jl","text":"using StatsBase\nusing PrettyTables\nusing Random\nusing CorrelationFunctions\nusing CorrelationTrackers\n\n# Create our system\nsystem = rand(MersenneTwister(348), 0:1, (20, 30))\n\n# Suppose we want to track $S_2^1(r)$ and $L_2^0(r)$.\ntracking = [TrackedData(Directional.s2, 1), \n            TrackedData(Directional.l2, 1)]\n\n# Create the tracker. It may require some time if your system is big.\ntracker = CorrelationTracker{Int, 2}(system; tracking = tracking)\n\n# Do some phase swapping in both the original system and the tracker\nsystem[12, 25]  = 1 - system[12, 25]\nsystem[11, 20]  = 1 - system[11, 20]\ntracker[12, 25] = 1 - tracker[12, 25]\ntracker[11, 20] = 1 - tracker[11, 20]\n\npretty_table(stdout, hcat(Directional.s2(system,  1) |> mean,\n                          Directional.l2(system,  1) |> mean,\n                          Directional.s2(tracker, 1) |> mean,\n                          Directional.l2(tracker, 1) |> mean);\n                     header = [\"S2 recalculated\", \"L2 recalculated\",\n                               \"S2 tracked\",      \"L2 tracked\"])","category":"page"},{"location":"index.html#API","page":"CorrelationTrackers.jl","title":"API","text":"","category":"section"},{"location":"index.html","page":"CorrelationTrackers.jl","title":"CorrelationTrackers.jl","text":"TrackedData\ndefault_trackers\nCorrelationTracker\ntracked_data\nsoftupdate","category":"page"},{"location":"index.html#CorrelationTrackers.TrackedData","page":"CorrelationTrackers.jl","title":"CorrelationTrackers.TrackedData","text":"TrackedData{T}(func :: Function, phase :: T)\n\nConstruct a pair of correlation function func and phase phase which must be tracked in CorrelationTracker.\n\nExamples\n\njulia> TrackedData(Directional.l2, 1)\nTrackedData{Int64}(CorrelationFunctions.Directional.l2, 1)\n\n\n\n\n\n","category":"type"},{"location":"index.html#CorrelationTrackers.default_trackers","page":"CorrelationTrackers.jl","title":"CorrelationTrackers.default_trackers","text":"default_trackers(T)\n\nConstruct a vector of correlation functions which are tracked by default (that is S_2^1(x), L_2^1(x) and L_2^0(x)). T is the type of x.\n\n\n\n\n\n","category":"function"},{"location":"index.html#CorrelationTrackers.CorrelationTracker","page":"CorrelationTrackers.jl","title":"CorrelationTrackers.CorrelationTracker","text":"CorrelationTracker{T, N}(system   :: AbstractArray{T, N}; \n                         tracking = default_trackers(T),\n                         periodic = false[, directions][, kwargs...])\n\nCreate correlation functions tracker.\n\nCreate correlation tracker for the array system. tracking is a vector of TrackedData structures which specify correlation functions you wish to track. periodic and direction have the same meaning as in the most functions in CorrelationFunctions.jl package. Additional arguments such as len may be passed in kwargs.\n\nReturned tracker supports interface of AbstractArray (e.g. you can perform element-wise read and write operations).\n\nExamples\n\njulia> let\n       system = rand(MersenneTwister(35), 0:1, (30, 10))\n       tracker = CorrelationTracker{Int,2}(system)\n       end\n30×10 CorrelationTracker{Int64, 2}:\n 0  1  0  1  1  0  0  1  1  0\n 1  1  1  0  0  0  0  0  1  1\n 0  0  0  0  0  0  1  1  0  1\n 1  1  1  0  1  1  1  0  1  0\n 0  1  0  0  1  0  0  1  1  1\n 0  0  0  0  0  0  1  0  1  1\n 0  0  1  0  1  1  0  1  0  1\n 1  0  0  1  0  0  1  0  1  0\n 0  1  1  0  0  1  1  1  1  1\n 0  0  1  1  1  1  0  0  0  0\n 0  0  1  1  0  0  1  1  1  0\n 0  1  0  0  0  1  0  0  1  0\n 1  0  0  1  0  0  1  1  0  1\n 0  1  0  1  0  0  1  1  1  0\n 1  1  0  1  1  1  0  1  0  1\n 1  1  1  0  0  0  0  1  0  1\n 1  0  0  1  0  0  1  1  1  0\n 0  0  0  1  0  0  0  1  1  0\n 1  0  1  0  1  0  0  0  1  0\n 1  0  0  1  0  0  0  0  0  1\n 1  1  1  0  1  0  1  0  1  1\n 0  1  0  1  1  0  0  1  0  1\n 0  0  0  1  0  0  1  1  1  1\n 0  0  1  1  1  1  0  1  1  0\n 1  0  1  1  0  0  0  0  0  1\n 1  1  0  1  0  1  1  0  1  0\n 0  1  1  0  0  1  1  0  1  0\n 0  1  0  0  1  0  0  1  0  0\n 1  1  0  0  1  1  0  0  0  1\n 0  0  1  1  0  1  1  1  1  0\n\n\n\n\n\n","category":"type"},{"location":"index.html#CorrelationTrackers.tracked_data","page":"CorrelationTrackers.jl","title":"CorrelationTrackers.tracked_data","text":"tracked_data(x :: CorrelationTracker)\n\nReturn an iterator over TrackedData objects which are tracked by the tracker.\n\n\n\n\n\n","category":"function"},{"location":"index.html#CorrelationTrackers.softupdate","page":"CorrelationTrackers.jl","title":"CorrelationTrackers.softupdate","text":"softupdate(x :: CorrelationTracker, val, idx...)\n\nPerform soft update of correlation functions.\n\nThis function is equivalent to x[idx...] = val, but instead of modifying x it returns another CorrelationTracker with updated correlation functions. The new tracker is read-only, i.e you cannot assign values to elements of an underlying array.\n\nExamples\n\njulia> let\n       array = rand(0:1, (200, 200))\n       tracker = CorrelationTracker{Int, 2}(array)\n       su = softupdate(tracker, 1 - tracker[43, 102], 43, 102)\n       tracker[43, 102] = 1 - tracker[43, 102]\n       Directional.s2(tracker, 1)[:x] == Directional.s2(su, 1)[:x]\n       end\ntrue\n\n\n\n\n\n","category":"function"},{"location":"index.html#Caveats","page":"CorrelationTrackers.jl","title":"Caveats","text":"","category":"section"},{"location":"index.html","page":"CorrelationTrackers.jl","title":"CorrelationTrackers.jl","text":"Currently Directional.l2 and Directional.s2 return internal structures of supplied tracker. Do not modify them. This behavior can be changed to a safer one in the future.","category":"page"}]
}
