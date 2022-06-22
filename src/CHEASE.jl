"""
Julia wrapper for CHEASE fixed boundary equilibrium solver (https://gitlab.epfl.ch/spc/chease.git)
"""

module CHEASE

__precompile__(true)

using Fortran90Namelists
using EFIT
using Equilibrium
import Dates

μ_0 = 1.25663706212e-6

function run_chease(ϵ::Real, z_axis::Real, pressure_sep::Real, Bt_center::Real, r_center::Real, Ip::Real, r_bound::AbstractVector{<:Real}, z_bound::AbstractVector{<:Real}, mode::Int64, rho_psi::Union{Missing,AbstractVector{<:Real}}, pressure::AbstractVector{<:Real}, j_tor::AbstractVector{<:Real})

    # File path and directory creation
    chease_dir = joinpath(dirname(abspath(@__FILE__)), "..")
    template_dir = joinpath(chease_dir, "templates")

    executable = joinpath(chease_dir, "executables", "chease_m1_ARM_gfortran")
    chease_namelist = joinpath(template_dir, "chease_namelist_OMFIT")
    run_dir = mkdir(joinpath(chease_dir, "rundir", string(Dates.now())))

    cp(chease_namelist, joinpath(run_dir, "chease_namelist"))
    cd(run_dir)

    # Edit chease namelist
    edit_namelist(chease_namelist, Bt_center, r_center, Ip)

    # Create EQOUT file
    write_EXPEQ_file(ϵ, z_axis, pressure_sep, r_bound, z_bound, mode, rho_psi, pressure, j_tor)

    # run chease
    write("chease.output", read(`$(executable)`))

    # read output
    EFITEquilibrium = read_chease_output(joinpath(run_dir, "EQDSK_COCOS_01.OUT"))

    # return arrays
    return EFITEquilibrium
    
end

function edit_namelist(chease_namelist, Bt_center, r_center, Ip)
    nml = readnml(chease_namelist)
    eqdata = nml[:EQDATA]

    eqdata[:R0EXP] = abs(r_center)
    eqdata[:B0EXP] = abs(Bt_center)
    eqdata[:CURRT] = abs(Ip / (r_center * Bt_center / μ_0))
    eqdata[:SIGNB0XP] = sign(Bt_center)
    eqdata[:SIGNIPXP] = sign(Ip)

    writenml(joinpath(pwd(), "chease_namelist"), nml; verbose=false)

    # Annoyingly CHEASE cares about the order of the namelist items...
    open(joinpath(pwd(), "chease_namelist")) do f
        encounter = false
        list_1 = [""]
        list_2 = [""]
        while ! eof(f)
            line = readline(f)
            if !encounter && line !== "/"
    
                list_1 = hcat(list_1,"$line")
            elseif encounter
                list_2 = hcat(list_2,"$line")
            end
    
            if line == "/"
                encounter = true
            end
        end
        open("chease_namelist", "w") do file
            for line in hcat(list_2,list_1, "/")
                write(file, line, "\n")
            end
        end
    end

end

function write_EXPEQ_file(ϵ, z_axis, pressure_sep, r_bound, z_bound, mode, rho_psi, pressure, j_tor)
    write_list = [string(ϵ), string(z_axis), string(pressure_sep)]
    @assert length(r_bound) == length(z_bound) "R,Z boundary arrays must have the same shape"
    write_list = vcat(write_list, string(length(r_bound)))
    for (r,z) in zip(r_bound, z_bound)
        write_list = vcat(write_list, "$r    $z")
    end
    @assert length(rho_psi) == length(pressure) == length(j_tor) "rho_psi, presssure and j_tor arrays must have the same shape"
    write_list = vcat(write_list, "$(length(pressure))    $(string(mode)[1])")
    write_list = vcat(write_list, "$(string(mode)[2])    0")
    write_list = vcat(write_list, map(string, rho_psi))
    write_list = vcat(write_list, map(string, pressure))
    write_list = vcat(write_list, map(string, j_tor))

    touch("EXPEQ")
    open("EXPEQ" , "w") do file
        for line in write_list
            write(file, "$line \n")
        end
    end
end


function read_chease_output(EQDSK)
    return efit(readg(EQDSK), 1)
end

export run_chease

end # module