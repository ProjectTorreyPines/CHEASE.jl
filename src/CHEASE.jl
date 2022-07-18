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
    function run_chease(
        ϵ::Real,
        z_axis::Real,
        pressure_sep::Real,
        Bt_center::Real,
        r_geo::Real,
        Ip::Real,
        r_bound::AbstractVector{<:Real},
        z_bound::AbstractVector{<:Real},
        mode::Int64,
        rho_psi::Union{Missing,AbstractVector{<:Real}},
        pressure::AbstractVector{<:Real},
        j_tor::AbstractVector{<:Real})

This function executes chease given the above set-up and handles the file-io
Returns an EFITEquilibrium struct (see Equilibrium/src/efit.jl)
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
    clear_workdir::Bool)
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
    cd(run_dir)

    # Edit chease namelist
    edit_chease_namelist(chease_namelist, Bt_center, r_geo, Ip, r_bound[1:end-1], z_bound[1:end-1])

    # Create EQOUT file
    write_EXPEQ_file(ϵ, z_axis, pressure_sep, r_geo, Bt_center, Ip, r_bound[1:end-1], z_bound[1:end-1], mode, rho_psi, pressure, j_tor)

    # run chease
    write("chease.output", read(`$(executable)`))

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