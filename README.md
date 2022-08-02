# CHEASE.jl
### Julia wrapper for the CHEASE fixed boundary equilibrium solver (https://gitlab.epfl.ch/spc/chease.git)

see `examples/run_chease_example.jl` on how to use this package.

Exectues chease using the following scalars and arrays:
```julia
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
        clear_workdir::Bool=true,
        extra_box_fraction::Real=0.33)
```
