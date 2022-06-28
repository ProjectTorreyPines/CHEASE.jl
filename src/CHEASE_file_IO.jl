μ_0 = 1.25663706212e-6

"""
    write_EXPEQ_file(ϵ, z_axis, pressure_sep, r_bound, z_bound, mode, rho_psi, pressure, j_tor)

This function writes a EXPEQ file for CHEASE given the above arrays and scalars
"""
function write_EXPEQ_file(ϵ, z_axis, pressure_sep, r_center, Bt_center, r_bound, z_bound, mode, rho_psi, pressure, j_tor)

    # Normalize from SI to chease units
    pressure_sep_norm = pressure_sep / (Bt_center^2 / μ_0)
    pressure_norm = pressure / (Bt_center^2 / μ_0)
    j_tor_norm = j_tor / (Bt_center / (r_center * μ_0))

    r_bound_norm = r_bound / r_center
    z_bound_norm = z_bound / r_center

    write_list = [string(ϵ), string(z_axis), string(pressure_sep_norm)]
    @assert length(r_bound) == length(z_bound) "R,Z boundary arrays must have the same shape"
    write_list = vcat(write_list, string(length(r_bound)))
    for (r,z) in zip(r_bound_norm, z_bound_norm)
        write_list = vcat(write_list, "$r    $z")
    end
    @assert length(rho_psi) == length(pressure) == length(j_tor) "rho_psi, presssure and j_tor arrays must have the same shape"
    write_list = vcat(write_list, "$(length(pressure))    $(string(mode)[1])")
    write_list = vcat(write_list, "$(string(mode)[2])    0")
    write_list = vcat(write_list, map(string, rho_psi))
    write_list = vcat(write_list, map(string, pressure_norm))
    write_list = vcat(write_list, map(string, j_tor_norm))

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
function edit_chease_namelist(chease_namelist, Bt_center, r_center, Ip, r_bound, z_bound)
    nml = readnml(chease_namelist)
    eqdata = nml[:EQDATA]

    eqdata[:R0EXP] = abs(r_center)
    eqdata[:B0EXP] = abs(Bt_center)
    eqdata[:CURRT] = abs(Ip / (r_center * Bt_center / μ_0))
    eqdata[:SIGNB0XP] = sign(Bt_center)
    eqdata[:SIGNIPXP] = sign(Ip)

#    eqdata[:COCOS_IN] = 11

    # box length
    r_max = maximum(r_bound)
    r_min = minimum(r_bound)
    r_mid = 0.5 * (r_max + r_max)
    z_max = maximum(z_bound)
    z_min = minimum(z_bound)

    eqdata[:RBOXLEN] = 1.25 * (r_max - r_min)
    eqdata[:RBOXLFT] = r_mid - 0.5 * eqdata[:RBOXLEN]
    eqdata[:ZBOXLEN] = 1.25 * (z_max - z_min)
    eqdata[:ZBOXMID] = 0.5 * (z_max + z_min)

    writenml(joinpath(pwd(), "chease_namelist"), nml; verbose=false)

end

"""
    read_chease_output(EQDSK)

This function reads the EQDSK output file from chease using the EFIT package and returns an EFITEquilibrium
"""
function read_chease_output(EQDSK)
    return readg(EQDSK)
end