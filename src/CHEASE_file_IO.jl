μ_0 = 1.25663706212e-6

"""
    write_EXPEQ_file(ϵ, z_axis, pressure_sep, r_center, Bt_center, Ip, r_bound, z_bound, mode, rho_pol, pressure, j_tor)

This function writes a EXPEQ file for CHEASE given the above arrays and scalars
"""
function write_EXPEQ_file(ϵ, z_axis, pressure_sep, r_center, Bt_center, Ip, r_bound, z_bound, mode, rho_pol, pressure, j_tor)

    # Normalize from SI to chease units
    pressure_sep_norm = pressure_sep / (Bt_center^2 / μ_0)
    pressure_norm = pressure / (Bt_center^2 / μ_0)
    j_tor_norm = abs.(j_tor / (Bt_center / (r_center * μ_0)))

    ip_sign = sign(Ip)
    bt_sign = sign(Bt_center)
    if ip_sign == 1 && bt_sign == 1
        j_tor_norm .*= 1
    elseif ip_sign == 1 && bt_sign == -1
        j_tor_norm .*= 1
    elseif ip_sign == -1 && bt_sign == -1
        j_tor_norm .*= 1
    else
        j_tor_norm .*= -1
    end

    r_bound_norm = r_bound / r_center
    z_bound_norm = z_bound / r_center

    write_list = [string(ϵ), string(z_axis), string(pressure_sep_norm)]
    @assert length(r_bound) == length(z_bound) "R,Z boundary arrays must have the same shape"
    write_list = vcat(write_list, string(length(r_bound)))
    for (r, z) in zip(r_bound_norm, z_bound_norm)
        write_list = vcat(write_list, "$r    $z")
    end
    @assert length(rho_pol) == length(pressure) == length(j_tor) "rho_pol, presssure and j_tor arrays must have the same shape"
    write_list = vcat(write_list, "$(length(pressure))    $(string(mode)[1])")
    write_list = vcat(write_list, "$(string(mode)[2])    0")
    write_list = vcat(write_list, map(string, rho_pol))
    write_list = vcat(write_list, map(string, pressure_norm))
    write_list = vcat(write_list, map(string, j_tor_norm))

    touch("EXPEQ")
    open("EXPEQ", "w") do file
        for line in write_list
            write(file, "$line \n")
        end
    end
end

"""
    write_chease_namelist(chease_namelist, Bt_center, r_center, Ip, r_bound, z_bound; extra_box_fraction=0.33)

This function writes the chease_namelist using the Fortran90Namelists package
"""
function write_chease_namelist(chease_namelist, Bt_center, r_center, Ip, r_bound, z_bound; rescale_eq_to_ip=false, extra_box_fraction=0.33)
    nml = readnml(chease_namelist)
    eqdata = nml[:EQDATA]

    eqdata[:R0EXP] = r_center
    eqdata[:B0EXP] = Bt_center
    eqdata[:CURRT] = abs(Ip / (r_center * Bt_center / μ_0))
    eqdata[:SIGNB0XP] = sign(Bt_center)
    eqdata[:SIGNIPXP] = sign(Ip)
    eqdata[:NT] = 80
    eqdata[:COCOS_IN] = 11
    if rescale_eq_to_ip
        eqdata[:NCSCAL] = 2
    else
        eqdata[:NCSCAL] = 4
    end
    eqdata[:NPROPT] = -2
    eqdata[:NPPFUN] = 8
    eqdata[:EPSLON] = 1e-9

    # box length
    r_max = maximum(r_bound)
    r_min = minimum(r_bound)
    z_max = maximum(z_bound)
    z_min = minimum(z_bound)

    r_extra = (r_max - r_min) * extra_box_fraction
    z_extra = (z_max - z_min) * extra_box_fraction
    eqdata[:RBOXLEN] = (r_max - r_min) + r_extra * 2.0
    eqdata[:RBOXLFT] = max(0, r_min - r_extra)
    eqdata[:ZBOXLEN] = (z_max - z_min) + z_extra * 2.0
    eqdata[:ZBOXMID] = (z_max + z_min) / 2.0

    writenml(joinpath(pwd(), "chease_namelist"), nml; verbose=false)

end

"""
    read_chease_output(EQDSK)

This function reads the EQDSK output file from chease using the EFIT package and returns an EFITEquilibrium
"""
function read_chease_output(EQDSK)
    return MXHEquilibrium.readg(EQDSK)
end