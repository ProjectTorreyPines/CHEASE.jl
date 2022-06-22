μ_0 = 1.25663706212e-6

"""
    write_EXPEQ_file(ϵ, z_axis, pressure_sep, r_bound, z_bound, mode, rho_psi, pressure, j_tor)

This function writes a EXPEQ file for CHEASE given the above arrays and scalars
"""
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

"""
    edit_chease_namelist(chease_namelist, Bt_center, r_center, Ip)

This function edits the chease_namelist using the Fortran90Namelists package
"""
function edit_chease_namelist(chease_namelist, Bt_center, r_center, Ip)
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

"""
    read_chease_output(EQDSK)

This function reads the EQDSK output file from chease using the EFIT package and returns an EFITEquilibrium
"""
function read_chease_output(EQDSK)
    return efit(readg(EQDSK), 1)
end