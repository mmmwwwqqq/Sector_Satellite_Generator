include("materials.jl")
include("utils.jl")

mutable struct Gearing
    z_p::Int32
    a_p::Float64
    e::Float64
    d_p::Float64
    s::Int8
    z_c::Int32
    λ::Float64
    satellite_material::Union{Nothing, Material}
    pin_material::Union{Nothing, Material}

    function Gearing(;z_p::Int, a_p::Real, e::Real, d_p::Real, s::Int)
        λ = 2 * z_p * e / a_p
        z_c = z_p - s
        return new(z_p, a_p, e, d_p, s, z_c, λ, nothing, nothing)
    end

    function GearingByLambda(;z_p, a_p, λ, d_p, s)
        e = λ * a_p / (2 * z_p)
        z_c = z_p - s 
        return new(z_p, a_p, e, d_p, s, z_c, λ, nothing, nothing)
    end
end

function cycloid_point(t, gear)
    a_p = gear.a_p
    e = gear.e
    z_p = gear.z_p
    s = gear.s
    Cx = a_p / 2 * sin(t) - e * sin(z_p * t)
    Cy = a_p / 2 * cos(t) - e * s * cos(z_p * t)
    return [Cx, Cy]
end

# function cycloid_points(gear; N=10000)
#     for i in 1:N


function normal(t, gear)
    λ = gear.λ
    s = gear.s
    z_p = gear.z_p
    sqr = sqrt(1 + λ^2 - 2 * λ * cos((z_p - s)*t))
    Nx = s * (-sin(t) + λ * s * sin(z_p * t)) / sqr
    Ny = s * (-cos(t) + λ * cos(z_p * t)) / sqr
    return [Nx, Ny]
end

function satellite_point(t, gear)
    Cx, Cy = cycloid_point(t, gear)
    Nx, Ny = normal(t, gear)
    d_p = gear.d_p
    Px = Cx + d_p / 2 * Nx
    Py = Cy + d_p / 2 * Ny
    return [Px, Py]
end

function satellite_border_rectangle(center_angle, satellite_angle, satellite_height, gear)
    right_bottom = satellite_point(center_angle + satellite_angle / 2, gear)
    left_bottom = satellite_point(center_angle - satellite_angle / 2, gear)
    dx = satellite_height * sin(center_angle)
    dy = satellite_height * cos(center_angle)
    left_upper = left_bottom .+ [dx, dy]
    right_upper = right_bottom .+ [dx, dy]
    border_coordinates = zeros(4, 2)
    for i in 1:2
        border_coordinates[:, i] .= [right_bottom[i], right_upper[i], left_upper[i], left_bottom[i]]
    end
    return border_coordinates
end

function pin_center(t, number, gear)

    t_center = 2π / z_p * number + t
    center_x, center_y = cycloid_point(t_center, gear)

    return [center_x, center_y]
end

function pin_array(t, number, gear)
    center_x, center_y = pin_center(t, number, gear)
    t_circle = range(0, 2π, 100)
    x_circle = center_x .+ cos.(t_circle) * d_p / 2
    y_circle = center_y .+ sin.(t_circle) * d_p / 2
    return x_circle; y_circle
end

# function satellite

