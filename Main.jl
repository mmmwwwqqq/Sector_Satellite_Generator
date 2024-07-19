using Plots
using LinearAlgebra
using DelimitedFiles
using Gtk
include("gearing.jl")
include("utils.jl")
include("materials.jl")

# Конвертировать градусы в радианы
function deg2rad(deg)
    return deg * π / 180.0
end

# Функция для обновления графика и вывода значения λ
function update_plot(z_p, a_p, d_p, e, θ, satellite_center_angle, β, z_s, N, number_satellite_points, s, xlim_range, ylim_range)
    cycloid_points = zeros(N, 2)
    satellite_points = zeros(number_satellite_points, 2)

    gear = Gearing(z_p=z_p, a_p=a_p, e=e, d_p=d_p, s=s)
    gear.pin_material = mSteel40
    gear.satellite_material = mKaprolonPA16
    lambda_val = gear.λ  # Получаем значение λ
    println("gearing lambda = $lambda_val")
    angle_between_excenters = 2π / z_s
    p = plot()

    # Выводим цевки
    for i in 1:z_p
        t_center = 2π * i / z_p - gear.s * β / z_p
        center_x, center_y = gear.a_p / 2 .* [sin(t_center), cos(t_center)]
        t_circle = range(0, 2π, length=N)
        x_circle = center_x .+ cos.(t_circle) * d_p / 2
        y_circle = center_y .+ sin.(t_circle) * d_p / 2
        plot!(p, x_circle, y_circle, color=:black, legend=false, aspect_ratio=:1)
    end

    all_contact_points = zeros(z_s, z_p, 2)
    all_normals = zeros(z_s, z_p, 2)
    all_masks = zeros(Bool, z_s, z_p)

    for (i, tau) in enumerate(range(0, 2π - angle_between_excenters, length=z_s))
        t_new_satellite = range(satellite_center_angle - θ / 2 + tau / z_p, satellite_center_angle + θ / 2 + tau / z_p, length=number_satellite_points)
        t_all_contact_points = zeros(z_p)
        for j in 1:z_p
            t_all_contact_points[j] = 2π * (j - 1) / z_p + tau / z_p - gear.s * β / z_p
        end
        anonim = x -> angle_is_between(x, satellite_center_angle - θ / 2 + tau / z_p, satellite_center_angle + θ / 2 + tau / z_p)
        mask_contact_points = anonim.(t_all_contact_points)
        if θ ≈ 2π || θ >= 2π
            mask_contact_points .= true
        end
        all_masks[i, :] .= mask_contact_points
        t_contact_points = t_all_contact_points[mask_contact_points]
        z_c = length(t_contact_points)
        A_01 = rotate_matrix_2d(gear.s * (tau)) * shift_matrix_2d([0; gear.s * e]) * rotate_matrix_2d(-gear.s * (tau))
        A_02 = rotate_matrix_2d(-tau / z_p)
        A_03 = rotate_matrix_2d(gear.s * (tau - gear.s * β)) * shift_matrix_2d([0; gear.s * e]) * rotate_matrix_2d(-gear.s * (tau - gear.s * β)) * inv(A_01)
        gear_center = [-gear.e * sin(tau + β); -gear.s * gear.e * cos(tau + β); 1]
        satellite_center = [0; 0; 1]
        satellite_center = A_02 * A_03 * A_01 * satellite_center
        gear_center = A_02 * A_03 * A_01 * gear_center
        satellite_coord_system_x = [0.5; 0; 1]
        satellite_coord_system_y = [0; 0.5; 1]
        satellite_coord_system_x = A_02 * A_03 * A_01 * satellite_coord_system_x
        satellite_coord_system_y = A_02 * A_03 * A_01 * satellite_coord_system_y

        for j in 1:z_p
            P = satellite_point(t_all_contact_points[j], gear)
            P = [P[1]; P[2]; 1]
            P = A_02 * A_03 * A_01 * P
            n = normal(t_all_contact_points[j], gear)
            n = [n[1]; n[2]; 1]
            n = A_02 * n
            all_contact_points[i, j, :] = P[1:2]
            all_normals[i, j, :] = n[1:2]
        end

        contact_points = zeros(z_c, 2)
        for j in 1:z_c
            cp = satellite_point(t_contact_points[j], gear)
            cp = [cp[1]; cp[2]; 1]
            cp = A_02 * A_03 * A_01 * cp
            contact_points[j, :] = cp[1:2]
        end

        plot!(p, contact_points[:, 1], contact_points[:, 2], label="Contact Points", xlabel="X", ylabel="Y", title="Contact Points Scatter Plot", seriestype=:scatter)

        new_satellite_points = zeros(number_satellite_points, 2)
        out_satellite_points = zeros(number_satellite_points, 2*z_s)
        for j in 1:number_satellite_points
            new_satellite_points[j, :] .= satellite_point(t_new_satellite[j], gear)
            P = [new_satellite_points[j, 1]; new_satellite_points[j, 2]; 1]
            P = A_02 * A_03 * A_01 * P
            new_satellite_points[j, :] .= P[1:2]
        end

        plot!(p, new_satellite_points[:, 1], new_satellite_points[:, 2], label="Satellite Points", xlabel="X", ylabel="Y", title="Satellite Points Scatter Plot")
        
        # Debugging statements to check file writing
        try
            output_file_path = joinpath(@__DIR__, "output_$(i).txt")
            println("Saving file to $output_file_path")
            writedlm(output_file_path, new_satellite_points, '\t')
            println("File saved successfully to $output_file_path")
        catch e
            println("Failed to save file: $e")
        end
    end

    plot!()

    # Устанавливаем новые границы графика
    xlims!(p, (xlim_range[1], xlim_range[2]))
    ylims!(p, (ylim_range[1], ylim_range[2]))

    Plots.savefig(p, "Plot_Gear1.svg")
    return lambda_val  # Возвращаем значение lambda для вывода в интерфейсе
end

function save_points()
    try
        output_file_path = joinpath(@__DIR__, "output_$(i).txt")
        println("Saving file to $output_file_path")
        writedlm(output_file_path, new_satellite_points, '\t')
        println("File saved successfully to $output_file_path")
    catch e
        println("Failed to save file: $e")
    end
end

# Функция для создания окна интерфейса
function create_window()
    win = GtkWindow("Генератор сателлитов", 800, 600)
    
    hbox = Gtk.Box(:h, 5)
    push!(win, hbox)
    
    vbox = Gtk.Box(:v, 5)
    push!(hbox, vbox)
    
    label_zp = Gtk.Label("z_p:")
    entry_zp = Gtk.Entry()
    set_gtk_property!(entry_zp, :text, "80")
    push!(vbox, label_zp)
    push!(vbox, entry_zp)

    label_ap = Gtk.Label("a_p:")
    entry_ap = Gtk.Entry()
    set_gtk_property!(entry_ap, :text, "400.0")
    push!(vbox, label_ap)
    push!(vbox, entry_ap)
    
    label_dp = Gtk.Label("d_p:")
    entry_dp = Gtk.Entry()
    set_gtk_property!(entry_dp, :text, "10")
    push!(vbox, label_dp)
    push!(vbox, entry_dp)

    label_e = Gtk.Label("e:")
    entry_e = Gtk.Entry()
    set_gtk_property!(entry_e, :text, "1.8")
    push!(vbox, label_e)
    push!(vbox, entry_e)

    label_θ = Gtk.Label("θ (degrees):")
    entry_θ = Gtk.Entry()
    set_gtk_property!(entry_θ, :text, "30.0")  # 30 градусов
    push!(vbox, label_θ)
    push!(vbox, entry_θ)

    label_satellite_center_angle = Gtk.Label("Satellite_center_angle (degrees):")
    entry_satellite_center_angle = Gtk.Entry()
    set_gtk_property!(entry_satellite_center_angle, :text, "45.0")  # 45 градусов
    push!(vbox, label_satellite_center_angle)
    push!(vbox, entry_satellite_center_angle)

    label_β = Gtk.Label("β (degrees):")
    entry_β = Gtk.Entry()
    set_gtk_property!(entry_β, :text, "0.0")
    push!(vbox, label_β)
    push!(vbox, entry_β)

    label_zs = Gtk.Label("z_s:")
    entry_zs = Gtk.Entry()
    set_gtk_property!(entry_zs, :text, "4")
    push!(vbox, label_zs)
    push!(vbox, entry_zs)
    
    label_N = Gtk.Label("N:")
    entry_N = Gtk.Entry()
    set_gtk_property!(entry_N, :text, "100")
    push!(vbox, label_N)
    push!(vbox, entry_N)

    label_number_satellite_points = Gtk.Label("number_satellite_points:")
    entry_number_satellite_points = Gtk.Entry()
    set_gtk_property!(entry_number_satellite_points, :text, "1000")
    push!(vbox, label_number_satellite_points)
    push!(vbox, entry_number_satellite_points) 

    label_s = Gtk.Label("s:")
    entry_s = Gtk.Entry()
    set_gtk_property!(entry_s, :text, "-1")
    push!(vbox, label_s)
    push!(vbox, entry_s)

    button_update = Gtk.Button("Update Plot")
    push!(vbox, button_update)

    # Добавляем кнопки для перемещения графики
    hbox_move = Gtk.Box(:h, 5)
    button_left = Gtk.Button("←")
    button_up = Gtk.Button("↑")
    button_down = Gtk.Button("↓")
    button_right = Gtk.Button("→")
    
    push!(hbox_move, button_left)
    push!(hbox_move, button_up)
    push!(hbox_move, button_down)
    push!(hbox_move, button_right)
    push!(vbox, hbox_move)
    
    image = Gtk.Image("Plot_Gear1.svg")
    push!(hbox, image)

    # Добавляем метку для вывода значения lambda
    label_lambda = Gtk.Label("Gearing Lambda:")
    push!(vbox, label_lambda)
    
    xlim_range = [-500, 500]
    ylim_range = [-500, 500]
    
    function on_button_clicked(widget)
        z_p = parse(Int, Gtk.get_gtk_property(entry_zp, :text, String))
        a_p = parse(Float64, Gtk.get_gtk_property(entry_ap, :text, String))
        d_p = parse(Int, Gtk.get_gtk_property(entry_dp, :text, String))
        e = parse(Float64, Gtk.get_gtk_property(entry_e, :text, String))
        θ = deg2rad(parse(Float64, Gtk.get_gtk_property(entry_θ, :text, String)))
        satellite_center_angle = deg2rad(parse(Float64, Gtk.get_gtk_property(entry_satellite_center_angle, :text, String)))
        β = deg2rad(parse(Float64, Gtk.get_gtk_property(entry_β, :text, String)))
        z_s = parse(Int, Gtk.get_gtk_property(entry_zs, :text, String))
        N = parse(Int, Gtk.get_gtk_property(entry_N, :text, String))
        number_satellite_points = parse(Int, Gtk.get_gtk_property(entry_number_satellite_points, :text, String))
        s = parse(Int, Gtk.get_gtk_property(entry_s, :text, String))
        
        lambda_val = update_plot(z_p, a_p, d_p, e, θ, satellite_center_angle, β, z_s, N, number_satellite_points, s, xlim_range, ylim_range)
        
        set_gtk_property!(label_lambda, :label, "Gearing Lambda: $(lambda_val)")  # Обновляем значение метки
        
        set_gtk_property!(image, :file, "Plot_Gear1.svg")
    end
    
    function on_button_move(widget, direction)
        step = 50  # Шаг перемещения
        if direction == "up"
            ylim_range[1] -= step
            ylim_range[2] -= step
        elseif direction == "down"
            ylim_range[1] += step
            ylim_range[2] += step
        elseif direction == "left"
            xlim_range[1] += step
            xlim_range[2] += step
        elseif direction == "right"
            xlim_range[1] -= step
            xlim_range[2] -= step
        end
        
        # Перестраиваем график с новыми границами
        on_button_clicked(widget)
    end

    function on_window_closed(widget)
        if isfile("Plot_Gear1.svg")
            rm("Plot_Gear1.svg")
        end
    end

    signal_connect(on_button_clicked, button_update, :clicked)
    signal_connect(button_left, :clicked) do widget on_button_move(widget, "left") end
    signal_connect(button_up, :clicked) do widget on_button_move(widget, "up") end
    signal_connect(button_down, :clicked) do widget on_button_move(widget, "down") end
    signal_connect(button_right, :clicked) do widget on_button_move(widget, "right") end
    signal_connect(on_window_closed, win, :destroy)

    showall(win)
end

println("ЗАПУСКАЕТСЯ!!!!!!!!!!!!!")
create_window()
