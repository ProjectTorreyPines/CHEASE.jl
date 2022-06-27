"""
Julia wrapper for CHEASE fixed boundary equilibrium solver (https://gitlab.epfl.ch/spc/chease.git)
"""

module CHEASE

__precompile__(true)

using Fortran90Namelists
import EFIT:readg
using Equilibrium
import Dates

# include CHEASE file handling functions
include("CHEASE_file_IO.jl")

"""
    function run_chease(
        ϵ::Real,
        z_axis::Real,
        pressure_sep::Real,
        Bt_center::Real,
        r_center::Real,
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
    r_center::Real,
    Ip::Real,
    r_bound::AbstractVector{<:Real},
    z_bound::AbstractVector{<:Real},
    mode::Int64,
    rho_psi::Union{Missing,AbstractVector{<:Real}},
    pressure::AbstractVector{<:Real},
    j_tor::AbstractVector{<:Real};
    keep_output=true)

    # File path and directory creation
    chease_dir = joinpath(dirname(abspath(@__FILE__)), "..")
    template_dir = joinpath(chease_dir, "templates")
    executable() = try
        return strip(read(`which chease`, String))
    catch
        return joinpath(chease_dir, "executables", "chease_m1_ARM_gfortran")
    end

    chease_namelist = joinpath(template_dir, "chease_namelist_OMFIT")
    run_dir = mkdir(joinpath(chease_dir, "rundir", string(Dates.now())))

    cp(chease_namelist, joinpath(run_dir, "chease_namelist"))
    cd(run_dir)

    # Edit chease namelist
    edit_chease_namelist(chease_namelist, Bt_center, r_center, Ip, r_bound, z_bound)

    # Create EQOUT file
    write_EXPEQ_file(ϵ, z_axis, pressure_sep, r_center, Bt_center, r_bound, z_bound, mode, rho_psi, pressure, j_tor)

    # run chease
    write("chease.output", read(`$(executable())`))

    # read output
    GEQDSKFile = read_chease_output(joinpath(run_dir, "EQDSK_COCOS_01.OUT"))

    if !keep_output
        rm(run_dir, force=true, recursive=true)
    end

    return GEQDSKFile
end

export run_chease

end # module