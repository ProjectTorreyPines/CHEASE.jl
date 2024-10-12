const μ_0 = 4pi * 1E-7

"""
    write_EXPEQ_file(
        ϵ::Float64,
        z_axis::Float64,
        pressure_sep::Float64,
        r_center::Float64,
        Bt_center::Float64,
        Ip::Float64,
        r_bound::Vector{Float64},
        z_bound::Vector{Float64},
        mode::Int,
        rho_pol::Vector{Float64},
        pressure::Vector{Float64},
        j_tor::Vector{Float64})

Writes a EXPEQ file for CHEASE given the above arrays and scalars
"""
function write_EXPEQ_file(
    ϵ::Float64,
    z_axis::Float64,
    pressure_sep::Float64,
    r_center::Float64,
    Bt_center::Float64,
    Ip::Float64,
    r_bound::Vector{Float64},
    z_bound::Vector{Float64},
    mode::Int,
    rho_pol::Vector{Float64},
    pressure::Vector{Float64},
    j_tor::Vector{Float64})

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

export write_EXPEQ_file
push!(document[:Base], :write_EXPEQ_file)

"""
    write_chease_namelist(
        Bt_center::Float64,
        r_center::Float64,
        Ip::Float64,
        r_bound::Vector{Float64},
        z_bound::Vector{Float64};
        rescale_eq_to_ip::Bool=false,
        extra_box_fraction::Float64=0.33)

Writes the chease namelist to the current folder
"""
function write_chease_namelist(
    Bt_center::Float64,
    r_center::Float64,
    Ip::Float64,
    r_bound::Vector{Float64},
    z_bound::Vector{Float64};
    rescale_eq_to_ip::Bool=false,
    extra_box_fraction::Float64=0.33)

    eqdata = Dict{Symbol,Any}()
    eqdata[:R0EXP] = r_center
    eqdata[:B0EXP] = Bt_center
    eqdata[:CURRT] = abs(Ip / (r_center * Bt_center / μ_0))
    eqdata[:SIGNB0XP] = sign(Bt_center)
    eqdata[:SIGNIPXP] = sign(Ip)
    eqdata[:NT] = 80 # number of theta points
    eqdata[:COCOS_IN] = 11
    if rescale_eq_to_ip
        eqdata[:NCSCAL] = 2
    else
        eqdata[:NCSCAL] = 4
    end
    eqdata[:NPROPT] = -2
    eqdata[:NPPFUN] = 8
    eqdata[:EPSLON] = 1e-6 # convergence
    eqdata[:RELAX] = 0.9

    # box size
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

    nml = """
&EQDATA
RELAX = $(eqdata[:RELAX])
NDIAGOP = 1
NBSEXPQ = 0
COCOS_IN = $(eqdata[:COCOS_IN])
COCOS_OUT = 1
NPROPT = $(eqdata[:NPROPT])
NIDEAL = 6
NPLOT = 1
NTCASE = 0
NSMOOTH = 1
NS = 40
NT = $(eqdata[:NT])
NPSI = 180
NCHI = 180
NISO = 180
NTNOVA = 12
CPRESS = 1.0
QSPEC = 0.7
CSSPEC = 0.0
CFNRESS = 1.0
NRSCAL = 0
NCSCAL = $(eqdata[:NCSCAL])
NTMF0 = 0
NBAL = 0
NBLOPT = 0
CFBAL = 10.0
NOPT = 0
R0EXP = $(eqdata[:R0EXP])
TENSPROF = -0.05
TENSBND = -0.05
NSURF = 6
ELONG = 2.045
TRIANG = 0.7
BEANS = 0.0
CETA = 0.24
SGMA = 0.0
ASPCT = 0.28
AT4 = 29500.0 -68768.0 272720.0 -1147400.0 2798300.0 -3873600.0 2842600.0 -852840.0
AT3 = 0.52503 0.92754 0.21896 -2.4078 8.1211 -13.87 11.653 -3.7942
AT2 = 1.5165 0.14189 -5.0417 36.759 -121.11 200.38 -162.23 51.152
ETAEI = 0.1
RPEOP = 0.5
RZION = 1.5
NPPFUN = $(eqdata[:NPPFUN])
NPP = 2
AP = 0.0 -0.8 2*0.0
NFUNC = 4
NSTTP = 2
AT = 0.0 -3.0761536 0.72318357 0.0
NSOUR = 8
NDIFPS = 0
NDIFT = 1
NMESHC = 1
NPOIDC = 2
SOLPDC = 0.7
CPLACE = 0.95 0.99 1.0
CWIDTH = 0.1 0.02 0.05
NMESHA = 0
NPOIDA = 1
SOLPDA = 0.1
APLACE = 0.0 0.7 1.0
AWIDTH = 0.05 0.07 0.05
NPOIDQ = 10
QPLACE = 2*1.0 2*2.0 2*3.0 2*4.0 2*4.41
QWIDTH = 0.13 0.04 0.09 0.04 0.07 0.02 0.04 2*0.01 0.001
NMESHD = 0
NPOIDD = 2
SOLPDD = 0.6
DPLACE = 2*-1.8 4.0
DWIDTH = 0.18 0.08 0.05
NMESHE = 0
NPOIDE = 4
SOLPDE = 0.5
EPLACE = 2*-1.7 2*1.7
EWIDTH = 0.18 0.08 0.18 0.08
EPSLON = $(eqdata[:EPSLON])
GAMMA = 1.6666666667
NTURN = 20
NBLC0 = 16
NPPR = 24
MSMAX = 1
NINMAP = 40
NINSCA = 40
NSYM = 0
NEGP = -1
NER = 1
NV = 40
NVEXP = 1
REXT = 10.0
R0W = 1.0
RZ0W = 0.0
NEQDXTPO = 1
NEQDSK = 0
PSISCL = 1.0
NRBOX = 129
NZBOX = 129
B0EXP = $(eqdata[:B0EXP])
CURRT = $(eqdata[:CURRT])
SIGNB0XP = $(eqdata[:SIGNB0XP])
SIGNIPXP = $(eqdata[:SIGNIPXP])
RBOXLEN = $(eqdata[:RBOXLEN])
RBOXLFT = $(eqdata[:RBOXLFT])
ZBOXLEN = $(eqdata[:ZBOXLEN])
ZBOXMID = $(eqdata[:ZBOXMID])
NIPR = 2
/
&NEWRUN
AL0 = -0.003
NWALL = 1
REXT = 10.0
NLGREN = .false.
WNTORE = 1.0
NV = 40
NLDIAG = 11*.true.
NAL0AUTO = 1
/
"""
    open(joinpath(pwd(), "chease_namelist"), "w") do file
        return write(file, nml)
    end

    return nml
end

export write_chease_namelist
push!(document[:Base], :write_chease_namelist)

"""
    read_chease_output(filename::String)

Reads gEQDSK output `filename` from CHEASE using the EFIT.jl package and returns an MXHEquilibrium object
"""
function read_chease_output(filename::String)
    return MXHEquilibrium.readg(filename; set_time=0.0)
end

export read_chease_output
push!(document[:Base], :read_chease_output)
