module CHEASE

import MXHEquilibrium
import EFIT

const document = Dict()
document[:Base] = Symbol[]
document[:IO] = Symbol[]

mutable struct Chease
    ϵ::Float64
    z_axis::Float64
    pressure_sep::Float64
    Bt_center::Float64
    r_geo::Float64
    Ip::Float64
    r_bound::Vector{Float64}
    z_bound::Vector{Float64}
    mode::Int
    rho_psi::Union{Missing,Vector{Float64}}
    pressure::Vector{Float64}
    j_tor::Vector{Float64}
    gfile::EFIT.GEQDSKFile
end

mutable struct MartianCHEASE
    ϵ::Float64
    z_axis::Float64
    pressure_sep::Float64
    Bt_center::Float64
    r_geo::Float64
    Ip::Float64
    r_bound::Vector{Float64}
    z_bound::Vector{Float64}
    mode::Int
    rho_psi::Union{Missing,Vector{Float64}}
    j_tor::Vector{Float64}
    pprime::Vector{Float64}
    number_walls::Int  # => NWBPS in EXPEQ
    wall_resistivity_type::Int # => NDATA in EXPEQ
    r_limiter::Union{Missing,Vector{Float64}}
    z_limiter::Union{Missing,Vector{Float64}}
end


# include CHEASE file handling functions
include("CHEASE_file_IO.jl")

"""
    run_chease(
        ϵ::Float64,
        z_axis::Float64,
        pressure_sep::Float64,
        Bt_center::Float64,
        r_geo::Float64,
        Ip::Float64,
        r_bound::Vector{Float64},
        z_bound::Vector{Float64},
        mode::Integer,
        rho_psi::Union{Missing,Vector{Float64}},
        pressure::Vector{Float64},
        j_tor::Vector{Float64};
        rescale_eq_to_ip::Bool=false,
        clear_workdir::Bool,
        extra_box_fraction::Float64=0.33)

This function executes chease given the above set-up and handles the file-io

Returns an EFITEquilibrium struct (see MXHEquilibrium/src/efit.jl)

The rescale_eq_to_ip option rescales the equilibrium to match Ip given (this is useful when using CHEASE from nothing where j_tor is madeup)
"""
function run_chease(
    ϵ::Float64,
    z_axis::Float64,
    pressure_sep::Float64,
    Bt_center::Float64,
    r_geo::Float64,
    Ip::Float64,
    r_bound::Vector{Float64},
    z_bound::Vector{Float64},
    mode::Integer,
    rho_psi::Union{Missing,Vector{Float64}},
    pressure::Vector{Float64},
    j_tor::Vector{Float64};
    rescale_eq_to_ip::Bool=false,
    clear_workdir::Bool,
    extra_box_fraction::Float64=0.33)

    # File path and directory creation
    chease_dir = joinpath(@__DIR__, "..")
    template_dir = joinpath(chease_dir, "templates")
    executable = try
        strip(read(`which chease`, String))
    catch
        if Sys.ARCH == :x86_64
            if Sys.islinux()
                joinpath(chease_dir, "executables", "chease_linux_x86")
            elseif Sys.isapple()
                joinpath(chease_dir, "executables", "chease_OSX_x86")
            else
                error("CHEASE.jl does not have a x86 binary for your OS")
            end
        elseif Sys.ARCH == :aarch64
            joinpath(chease_dir, "executables", "chease_m1_ARM")
        else
            error("CHEASE.jl does not have a binary for your CPU architecture")
        end
    end

    run_dir = mktempdir()
    if clear_workdir
        @debug "Running CHEASE in $run_dir"
    else
        @warn "CHEASE run directory $run_dir"
    end

    # NOTE: do not `cd` into run_dir. The process-wide working directory is
    # shared across all threads/tasks, so concurrent CHEASE runs would race
    # (capturing each other's run_dir as `old_dir` and failing once it is
    # cleaned up). Instead write every file with an absolute path under run_dir
    # and run the binary with its own working directory via `Cmd(...; dir=)`.
    output_path = joinpath(run_dir, "chease.output")
    try
        open(output_path, "w") do io
            # write chease namelist to file
            write_chease_namelist(Bt_center, r_geo, Ip, r_bound[1:end-1], z_bound[1:end-1]; rescale_eq_to_ip, extra_box_fraction, dir=run_dir)

            # Create EQOUT file
            write_EXPEQ_file(ϵ, z_axis, pressure_sep, r_geo, Bt_center, Ip, r_bound[1:end-1], z_bound[1:end-1], mode, rho_psi, pressure, j_tor; dir=run_dir)

            # run CHEASE
            return run(pipeline(Cmd(`$(executable)`; dir=run_dir); stdout=io, stderr=io))
        end

    catch e
        # show last 100 lines of chease.output file
        txt = open(output_path, "r") do io
            return split(read(io, String), "\n")
        end
        @error "ERROR running CHEASE\n...\n" * join(txt[max(1, length(txt) - 100):end], "\n")
        rethrow(e)
    end

    # read output
    gfile = read_chease_output(joinpath(run_dir, "EQDSK_COCOS_01.OUT"))

    if clear_workdir
        rm(run_dir; force=true, recursive=true)
    end

    # populate results data structure
    chease = Chease(
        ϵ,
        z_axis,
        pressure_sep,
        Bt_center,
        r_geo,
        Ip,
        r_bound,
        z_bound,
        mode,
        rho_psi,
        pressure,
        j_tor,
        gfile)

    return chease
end

export run_chease
push!(document[:Base], :run_chease)

end # module