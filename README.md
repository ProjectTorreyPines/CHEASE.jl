# CHEASE.jl
CHEASE implementation in Julia

see `examples/run_chease_example.jl on how to use this package.`

Exectues chease using the following scalars and arrays:
```julia
    function run_chease(
        ϵ::Real,
        z_axis::Real,
        pressure_sep::Real,
        Bt_center::Real,
        r_center::Real,
        Ip::Real,
        r_bound::AbstractVector{<:Real},
        z_bound::AbstractVector{<:Real},
        mode::Int64,
        rho_psi::Union{Missing,AbstractVector{<:Real}},
        pressure::AbstractVector{<:Real},
        j_tor::AbstractVector{<:Real})
```