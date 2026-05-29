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

function write_EXPEQ_file(eq::MartianCHEASE)

    pressure_sep_norm =
        eq.pressure_sep / (eq.Bt_center^2 / μ_0)

    pprime = 2 * pi * eq.pprime * eq.r_geo^2 * μ_0 / eq.Bt_center

    ## canNOT normalize here because MartianCHEASE uses FF', <Jtor> or <J//>
    #j_tor_norm =
    #    abs.(eq.j_tor ./ (eq.Bt_center/(eq.r_center*μ_0)))

    r_bound_norm = eq.r_bound ./ eq.r_geo
    z_bound_norm = eq.z_bound ./ eq.r_geo

    NWBPS = eq.number_walls
    NDATA = eq.wall_resistivity_type

    write_list = String[]

    push!(write_list,string(eq.ϵ))
    push!(write_list,string(eq.z_axis))
    push!(write_list,string(pressure_sep_norm))

    push!(
        write_list,
        string(length(eq.r_bound)," ",NWBPS, " ",NDATA)
    )

    for (r,z) in zip(r_bound_norm,z_bound_norm)
        push!(write_list,"$r    $z")
    end

    if NWBPS > 1 ## WHAT TO DO if > 2
        r_limiter_norm = eq.r_limiter ./ eq.r_geo
        z_limiter_norm = eq.z_limiter ./ eq.r_geo
        for (rw,zw) in zip(r_limiter_norm, z_limiter_norm)
            push!(write_list,"$(rw)    $(zw)")
        end
    
    end

    push!(write_list, "$(length(eq.pprime))")
    push!(write_list, "$(string(eq.mode))")
    
    append!(write_list,string.(eq.rho_psi))
    append!(write_list,string.(pprime))
    append!(write_list,string.(eq.j_tor))

    open("EXPEQ","w") do file
        for line in write_list
            write(file,"$line\n")
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

Base.@kwdef mutable struct CHEASEnamelist
    NEQDSK::Int     = 0
    NSURF::Int      = 6
    NTCASE::Int     = 0

    NBLOPT::Int     = 0
    NBSOPT::Int     = 0
    CPRESS::Float64 = 1.000
    CFBAL::Float64  = 1.0000 # set to 1. if NSCAL = 4

    NCSCAL::Int     = 2   # set to 4 if NOT scale q
    CSSPEC::Float64 = 0.000
    QSPEC::Float64  = 1.6185

    NTMF0::Int      = 0
    CURRT::Float64  = 0.3000

    NSTTP::Int      = 2
    NFUNC::Int      = 4
    NIPR::Int       = 1
    NISO::Int       = 100
    NIDEAL::Int     = 0

    NPPFUN::Int     = 4
    NPP::Int        = 1
    NPPR::Int       = 30

    NSOUR::Int      = 2
    NPROPT::Int     = 2

    NS::Int         = 60
    NT::Int         = 60
    NPSI::Int       = 240
    NCHI::Int       = 200

    NV::Int         = 160
    REXT::Float64   = 6.0
    NVEXP::Int      = 8
    R0W::Float64    = 0.90
    RZ0W::Float64   = 0.0

    NMESHA::Int     = 2
    SOLPDA::Float64 = 0.60
    QWIDTH0::Float64 = 0.30
    ROTE::Float64   = 0.0000
    NTOR::Int       = 1

    NPOIDQ::Int     = 6
    QSHAVE::Float64 = 100.0

    QPLACE::Vector{Float64} =
        [2.0000, 3.0000, 4.0000, 5.0000, 6.0000, 7.0000]

    QWIDTH::Vector{Float64} =
        [0.0009, 0.0007, 0.0006, 0.0006, 0.0005, 0.0005]

    NEGP::Int       = -1
    NER::Int        = 1

    EPSLON::Float64 = 1.0e-10
    GAMMA::Float64  = 1.6666666667

    MSMAX::Int      = 40
    NINMAP::Int     = 50
    NINSCA::Int     = 50

    NOPT::Int       = 0
    NPLOT::Int      = 1
    NBAL::Int       = 0

    B0EXP::Float64  = 1.5
    R0EXP::Float64  = 3.0
end

export CHEASEnamelist
push!(document[:Base], :CHEASEnamelist)

"""
    write_CHEASEnamelist(nl::CHEASEnamelist, filename::AbstractString="datain")

Writes a CHEASE namelist file from a CHEASEnamelist struct to `filename`.
"""
function write_CHEASEnamelist(
    nl::CHEASEnamelist,
    filename::AbstractString="datain"
)
    open(filename, "w") do io
        println(io, "***")
        println(io, "***    Example Torus")
        println(io, "***")
        println(io, "***")
        println(io, "&EQDATA")

        for field in fieldnames(CHEASEnamelist)
            val = getfield(nl, field)

            if val isa Vector
                println(io, "  $(field)(1) = ", join(val, ", "), ",")
            else
                println(io, "  $(field) = ", val, ",")
            end
        end

        println(io, "&END")
    end

    return filename
end

export write_CHEASEnamelist
push!(document[:Base], :write_CHEASEnamelist)
