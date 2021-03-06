using HDF5
using Distributed

struct FileWriter
    fid::HDF5.File
    channel::Union{Channel,RemoteChannel}
    total_trajectories::Int64
    timevector::Array{Float64,1}
    datasets::Dict
    writer::Task

    function FileWriter(filename::String, params::ModelParameters, total_trajectories::Integer, quantities)
        fid = h5open(filename, "w")

        timevector = collect(get_time(params))
        fid["t"] = timevector
        outpoints = length(timevector)

        paramdict = Dict(string(fn) => getfield(params, fn) for fn ∈ fieldnames(typeof(params)))

        for (par, val) in paramdict
            attributes(fid)[par] = val
        end

        datasets = Dict()
        # for quantity in ("FI", "QFI", "Jx", "Jy", "Jz", "Δjx", "Δjy", "Δjz", "xi2x", "xi2y", "xi2z", "FIstrong")
        for quantity in quantities
            # ds = d_create(fid, quantity, Float64, ((outpoints, total_trajectories), (outpoints, -1)), "chunk", (outpoints, 1))
            ds = create_dataset(fid, quantity, datatype(Float64), dataspace(outpoints, total_trajectories), chunk=(outpoints, 1))
            datasets[quantity] = ds
        end

        file_channel = RemoteChannel(() -> Channel{NamedTuple}(200))
        writer = @async writer_task(fid, datasets, file_channel)

        new(fid, file_channel, total_trajectories, timevector, datasets, writer)
    end
end

function writer_task(fid, datasets, channel)
    counter = 1
    attributes(fid)["stored_traj"] = 0
    @info "Writer ready!"
    while true
        try
            traj_result = take!(channel)
            for (d, data) in pairs(traj_result)
                if string(d) in keys(datasets)
                    try
                        ds = datasets[string(d)]
                        ds[:, counter] = data
                        flush(ds)
                    catch er
                        @error er
                    end
                end
            end
            @debug "Entry written to $(fid)"
            counter += 1
            write(attributes(fid)["stored_traj"], counter - 1)
            flush(fid)

        catch InvalidStateException
            @info "Channel closed"
            return fid
        end
    end
    return fid
end

function Base.close(fw::FileWriter)
    close(fw.channel)
    fid = fetch(fw.writer)
    close(fid)
    @info "FileWriter closed: $fid"
end

function Base.put!(fw::FileWriter, v)
    Base.put!(fw.channel, v)
end