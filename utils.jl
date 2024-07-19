function rotate_matrix_2d(phi)
    A = [cos(phi) sin(phi) 0;
        -sin(phi) cos(phi) 0;
        0 0 1;]
    return A
end

function shift_matrix_2d(px, py)
    A = [1 0 px;
         0 1 py;
         0 0 1]
    return A
end

function shift_matrix_2d(p)
    A = [1 0 p[1];
         0 1 p[2];
         0 0 1]
    return A
end

function clamp_angle(angle)
    new_angle = angle % 2π
    return (new_angle > 0) ? new_angle : 2π + new_angle
end

function clamp_angle!(angle)
    new_angle = angle % 2π
    angle = (new_angle > 0) ? new_angle : 2π + new_angle
end

function angle_is_between(angle, θ_left, θ_right)
    θ_1 = clamp_angle(θ_left)
    θ_2 = clamp_angle(θ_right)
    clamped_angle = clamp_angle(angle)
    if θ_1 > θ_2
        return clamped_angle >= θ_1 || clamped_angle <= θ_2
    end
    return clamped_angle >= θ_1 && clamped_angle <= θ_2
end

    