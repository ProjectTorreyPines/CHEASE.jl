# CHEASE.jl
### Julia wrapper for the CHEASE fixed boundary equilibrium solver (https://gitlab.epfl.ch/spc/chease.git)

see `examples/run_chease_example.jl` on how to use this package.

Exectues chease using the following scalars and arrays:
```julia
    function run_chease(
        ϵ::Real, # []
        z_axis::Real, #[m]
        pressure_sep::Real, #[Pa]
        Bt_center::Real,    #[T]
        r_center::Real,     #[m]
        Ip::Real,           #[A]
        r_bound::AbstractVector{<:Real}, #[m]
        z_bound::AbstractVector{<:Real}, #[m]
        mode::Int64, #[]
        rho_psi::Union{Missing,AbstractVector{<:Real}}, #[]
        pressure::AbstractVector{<:Real}, #[Pa] 
        j_tor::AbstractVector{<:Real})  #[A/m^2]
```