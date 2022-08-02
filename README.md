# CHEASE.jl
### Julia wrapper for the CHEASE fixed boundary equilibrium solver (https://gitlab.epfl.ch/spc/chease.git)

see `examples/run_chease_example.jl` on how to use this package.


## Call CHEASE in the following way:
(don't include the type information)

```julia
import CHEASE:run_chease

run_chease(
    ϵ::Real, # [-]
    z_axis::Real, # [m]
    pressure_sep::Real, # [Pa]
    Bt_center::Real, # [T]
    r_geo::Real, # [m]
    Ip::Real, # [A]
    r_bound::AbstractVector{<:Real}, # [m]
    z_bound::AbstractVector{<:Real}, # [m]
    mode::Integer, # [-]
    rho_psi::Union{Missing,AbstractVector{<:Real}}, # [-]
    pressure::AbstractVector{<:Real}, # [Pa]
    j_tor::AbstractVector{<:Real}; # [A.m^-2]
    rescale_eq_to_ip::Bool=false, 
    clear_workdir::Bool,
    extra_box_fraction::Real=0.33)
```
`rescale_eq_to_ip = true` runs CHEASE and rescales j_tor to match the Ip given

`rescale_eq_to_ip = false` runs CHEASE and tries to match pressure and j_tor
