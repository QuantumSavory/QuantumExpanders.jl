using TestItemRunner
using QuantumExpanders

# filter for the test
testfilter = ti -> begin
    exclude = Symbol[]

    if get(ENV, "JET_TEST", "") == "true"
        return :jet in ti.tags
    else
        push!(exclude, :jet)
    end

    if !(VERSION >= v"1.10")
        push!(exclude, :doctests)
        push!(exclude, :aqua)
    end

    if !(Base.Sys.islinux() & (Int===Int64))
        push!(exclude, :bitpack)
    end

    return all(!in(exclude), ti.tags)
end

println("Starting tests with $(Threads.nthreads()) threads out of `Sys.CPU_THREADS = $(Sys.CPU_THREADS)`...")

@run_package_tests filter=testfilter
