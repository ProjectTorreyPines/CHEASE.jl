"""
Julia wrapper for CHEASE fixed boundary equilibrium solver (https://gitlab.epfl.ch/spc/chease.git)
"""

module CHEASE

__precompile__(true)

using Fortran90Namelists
import Equilibrium
import EFIT

mutable struct Chease
    ϵ::Real
    z_axis::Real
    pressure_sep::Real
    Bt_center::Real
    r_geo::Real
    Ip::Real
    r_bound::AbstractVector{<:Real}
    z_bound::AbstractVector{<:Real}
    mode::Integer
    rho_psi::Union{Missing,AbstractVector{<:Real}}
    pressure::AbstractVector{<:Real}
    j_tor::AbstractVector{<:Real}
    gfile::EFIT.GEQDSKFile
end

# include CHEASE file handling functions
include("CHEASE_file_IO.jl")

"""
    run_chease(
        ϵ::Real,
        z_axis::Real,
        pressure_sep::Real,
        Bt_center::Real,
        r_geo::Real,
        Ip::Real,
        r_bound::AbstractVector{<:Real},
        z_bound::AbstractVector{<:Real},
        mode::Integer,
        rho_psi::Union{Missing,AbstractVector{<:Real}},
        pressure::AbstractVector{<:Real},
        j_tor::AbstractVector{<:Real};
        rescale_eq_to_ip::Bool=false, 
        clear_workdir::Bool,
        extra_box_fraction::Real=0.33)

This function executes chease given the above set-up and handles the file-io
Returns an EFITEquilibrium struct (see Equilibrium/src/efit.jl)
The rescale_eq_to_ip option rescales the equilibrium to match Ip given (This is useful when using CHEASE from nothing where j_tor is madeup)
"""
function run_chease(
    ϵ::Real,
    z_axis::Real,
    pressure_sep::Real,
    Bt_center::Real,
    r_geo::Real,
    Ip::Real,
    r_bound::AbstractVector{<:Real},
    z_bound::AbstractVector{<:Real},
    mode::Integer,
    rho_psi::Union{Missing,AbstractVector{<:Real}},
    pressure::AbstractVector{<:Real},
    j_tor::AbstractVector{<:Real};
    rescale_eq_to_ip::Bool=false,
    clear_workdir::Bool,
    extra_box_fraction::Real=0.33)

    # File path and directory creation
    chease_dir = joinpath(@__DIR__, "..")
    template_dir = joinpath(chease_dir, "templates")
    executable = try
        strip(read(`which chease`, String))
    catch
        joinpath(chease_dir, "executables", "chease_m1_ARM_gfortran")
    end

    chease_namelist = joinpath(template_dir, "chease_namelist_OMFIT")
    run_dir = mktempdir()
    @debug "Running CHEASE in $run_dir"

    cp(chease_namelist, joinpath(run_dir, "chease_namelist"))

    old_dir = pwd()
    try
        cd(run_dir)

        # Edit chease namelist
        write_chease_namelist(chease_namelist, Bt_center, r_geo, Ip, r_bound[1:end-1], z_bound[1:end-1]; rescale_eq_to_ip=rescale_eq_to_ip, extra_box_fraction=extra_box_fraction)

        # Create EQOUT file
        write_EXPEQ_file(ϵ, z_axis, pressure_sep, r_geo, Bt_center, Ip, r_bound[1:end-1], z_bound[1:end-1], mode, rho_psi, pressure, j_tor)

        # run chease
        open("chease.output", "w") do io
            run(pipeline(`$(executable)`; stdout=io, stderr=io))
        end

    catch
        # show last 100 lines of  chease.output
        txt = open("chease.output", "r") do io
            split(read(io, String), "\n")
        end
        @error "ERROR running CHEASE\n...\n" * join(txt[max(1, length(txt) - 100):end], "\n")

        rethrow()

    finally
        cd(old_dir)
    end

    # read output
    gfile = read_chease_output(joinpath(run_dir, "EQDSK_COCOS_01.OUT"))

    if clear_workdir
        rm(run_dir, force=true, recursive=true)
    else
        @warn "CHEASE run directory $run_dir"
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

end # module