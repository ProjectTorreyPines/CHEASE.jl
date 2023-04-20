"""
Julia wrapper for CHEASE fixed boundary equilibrium solver (https://gitlab.epfl.ch/spc/chease.git)
On OSX with M1 install gfortran https://github.com/R-macos/gcc-darwin-arm64/releases to use the CHEASE/executables/chease_m1_ARM executable
NOTE: was able to compile CHEASE with:
    gfortran -static-libgfortran -O3 -ffree-line-length-none -o chease itm_types.o euitm_schemas.o neobscoeffmod.o globals.o interpol.o prec_const.o sigmaneomod.o string_manipulation_tools.o checknanos_module.o cocos_module.o euitm_xml_parser.o assign_chease_codeparameters_reflist.o interpos_source.o chease_prog_effxml.o chease_effxml.o a_chease.o acopy.o aldlt.o apcoef.o apcoef2.o atcoef.o auxval.o away.o ballit.o baloon.o basis1.o basis2.o basis3.o basis4.o blines.o bltest.o bndspl.o bndfit.o bound.o bsexpeq.o bsfunc.o bstnzpro.o ccopy.o center.o check.o chipsi.o chipsimetrics.o cint.o conver.o copyap.o copyapp.o copyat.o cotrol.o cubrt.o curent.o cvzero.o direct.o drhodp.o dwy.o energy.o eqchease_mksa.o eqdim.o erdata.o errorch.o evlate.o fix_surface_near_axis.o four1.o fourfft.o fourier.o g_0.o g_1.o g_2.o g_3.o gauss.o gchi.o gdataext.o genout.o gijlin.o gloadd.o globals_init.o guess.o iarray.o identa.o identb.o indexx.o initia.o iodisk.o isamin.o ismax.o ismin.o isofind.o isofun.o isrchfge.o issum.o itipr.o ivar.o jnovaw.o labrun.o limita.o limitb.o ltxw.o lyv.o magaxe.o matrix.o mesage.o mesh.o metrictoitm_afterfpp.o msplcy.o mspline.o nerat.o nonlin.o norept.o ntridg.o oarray.o oldeq.o oldnew.o outgload.o outmksa.o outnvw.o outpen.o output.o outxt.o packme.o packmep.o page.o polyfun.o polynm.o ppbstr.o pprime.o pprm.o ppspln.o ppspln2.o premap.o preset.o prfunc.o priqqu.o prnorm.o profile.o psibox.o psicel.o psvol.o qplacs.o rarray.o realft.o reseti.o resetr.o resppr.o rmrad.o rscale.o runtim.o rvar.o rvar2.o rzbound.o scopyr.o setupa.o setupb.o shave.o smooth.o solovev.o solvit.o sort3.o splcy.o splcyp.o splifft.o spline.o ssum.o stchps.o stepon.o subsz.o surface.o surfadd.o surfrz.o surf_metrics_onaxis.o tcase.o test.o tetare.o tpsi.o tricyc.o tricycm.o tridagm.o tshift.o vacufft.o vacuum.o vacuumxt.o vlion.o vzero.o whtext.o witext.o wrtext.o wrtplot.o wrtmat.o wrtbin.o xtinit.o xtverifydims.o hamada.o neoart.o outgyro.o bscoeff.o outelit.o outastro.o ogyropsi.o prof2d_rz_to_fluxtheta.o interpos2d_cartesian.o write_ogyropsi.o gloqua.o mappin.o load_itm_dummy.o write_itm_dummy.o dpgbtrf_s.o

Instalation tip, install with latest gfortran
    Linux: conda install -c conda-forge gfortran
"""
module CHEASE

using Fortran90Namelists
import MXHEquilibrium
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
Returns an EFITEquilibrium struct (see MXHEquilibrium/src/efit.jl)
The rescale_eq_to_ip option rescales the equilibrium to match Ip given (this is useful when using CHEASE from nothing where j_tor is madeup)
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

    chease_namelist = joinpath(template_dir, "chease_namelist_OMFIT")
    run_dir = mktempdir()
    @debug "Running CHEASE in $run_dir"

    cp(chease_namelist, joinpath(run_dir, "chease_namelist"), force=true)

    old_dir = pwd()
    try
        cd(run_dir)

        # Edit chease namelist
        write_chease_namelist(chease_namelist, Bt_center, r_geo, Ip, r_bound[1:end-1], z_bound[1:end-1]; rescale_eq_to_ip, extra_box_fraction)

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