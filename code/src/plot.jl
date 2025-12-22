
function plot_unsteady_solid_body_rotation()
    # Figure 2a
    plot_convergence(
        [
            "../results/20250713_unsteady_solid_body_rotation_ec_N3/",
            "../results/20250713_unsteady_solid_body_rotation_es_N3/",
        ],
        plot_name = "unsteady_solid_body_rotation_convergence_N3.pdf",
        triangle_bottom_order = 4,
        triangle_top_order = 3,
    )
    # Figure 2b
    plot_convergence(
        [
            "../results/20250713_unsteady_solid_body_rotation_ec_N4/",
            "../results/20250713_unsteady_solid_body_rotation_es_N4/",
        ],
        plot_name = "unsteady_solid_body_rotation_convergence_N4.pdf",
        triangle_bottom_order = 5,
        triangle_top_order = 5,
    )

    # Figure 2c
    plot_convergence(
        [
            "../results/20250713_unsteady_solid_body_rotation_ec_p_refine/",
            "../results/20250714_unsteady_solid_body_rotation_es_p_refine/",
        ],
        triangle_bottom = false,
        triangle_top = false,
        plot_name = "unsteady_solid_body_rotation_p_refine_M4.pdf",
        xkey = "N",
        xlabel = L"Polynomial degree $N$",
        xticks = collect(2:10),
        xscale = identity,
        legend_position = (:left, :bottom),
    )
end

function plot_isolated_mountain()
    # Figure 3a
    plot_evolution(
        [
            "../results/20250713_well_balanced_ec/N3M20",
            "../results/20250713_well_balanced_es/N3M20",
        ],
        plot_name = "well_balanced_l2_h_evolution_N3M20.pdf",
        ykey = "l2_h",
        legend_position = (:left, :top),
        relative = false,
        ylabel = L"Normalized $L^2$ height error",
        xlims = [0, 15],
        ylims = [-0.2, 5],
        xticks = [0, 5, 10, 15],
        ynorm = 1e-14,
        exponent_text = L"\times 10^{-14}",
    )

    # Figure 3b
    plot_evolution(
        [
            "../results/20250713_isolated_mountain_ec/N3M20",
            "../results/20250713_isolated_mountain_es/N3M20",
        ],
        plot_name = "isolated_mountain_mass_evolution_N3M20.pdf",
        legend_position = (:left, :top),
        xlims = [0, 15],
        xticks = [0, 5, 10, 15],
        ylims = [-0.85, 4.5],
        ykey = "mass",
        ylabel = LaTeXString("Normalized mass change"),
        ynorm = 1e-14,
        exponent_text = L"\times 10^{-14}",
    )

    # Figure 3c
    plot_evolution(
        [
            "../results/20250713_isolated_mountain_ec/N3M20",
            "../results/20250713_isolated_mountain_es/N3M20",
        ],
        plot_name = "isolated_mountain_entropy_evolution_N3M20.pdf",
        legend_position = (:left, :bottom),
        xlims = [0, 15],
        xticks = [0, 5, 10, 15],
        ylims = [-8.25, 0.25],
        ykey = "entropy",
        ynorm = 1e-8,
        exponent_text = L"\times 10^{-8}",
    )
end

function plot_barotropic_instability()
    # Figure 5a
    plot_evolution(
        [
            "../results/20250525_steady_barotropic_instability_ec_M16/N3M16",
            "../results/20250525_steady_barotropic_instability_es_M16/N3M16",
        ],
        plot_name = "steady_barotropic_instability_l2_h_evolution_N3M16.pdf",
        ykey = "l2_h",
        legend_position = (:left, :top),
        relative = false,
        ynorm = 1e-3,
        exponent_text = L"\times 10^{-3}",
        ylabel = L"Normalized $L^2$ height error",
        xticks = [0, 3, 6, 9, 12],
        xlims = [0, 12],
        ylims = [-1, 12],
    )

    # Figure 5b
    plot_evolution(
        [
            "../results/20250520_steady_barotropic_instability_ec/N3M32",
            "../results/20250520_steady_barotropic_instability_es/N3M32",
        ],
        plot_name = "steady_barotropic_instability_l2_h_evolution_N3M32.pdf",
        ykey = "l2_h",
        legend_position = (:left, :top),
        relative = false,
        ynorm = 1e-3,
        exponent_text = L"\times 10^{-3}",
        ylabel = L"Normalized $L^2$ height error",
        xticks = [0, 3, 6, 9, 12],
        xlims = [0, 12],
        ylims = [-1, 12],
    )

    # Figure 5c
    plot_evolution(
        [
            "../results/20250520_steady_barotropic_instability_ec/N3M64",
            "../results/20250520_steady_barotropic_instability_es/N3M64",
        ],
        plot_name = "steady_barotropic_instability_l2_h_evolution_N3M64.pdf",
        ykey = "l2_h",
        legend_position = (:left, :top),
        relative = false,
        ynorm = 1e-3,
        exponent_text = L"\times 10^{-3}",
        ylabel = L"Normalized $L^2$ height error",
        xticks = [0, 3, 6, 9, 12],
        xlims = [0, 12],
        ylims = [-1, 12],
    )

    # Figure 5d
    plot_evolution(
        [
            "../results/20250525_steady_barotropic_instability_ec_M16/N3M16",
            "../results/20250525_steady_barotropic_instability_es_M16/N3M16",
        ],
        plot_name = "steady_barotropic_instability_entropy_evolution_N3M16.pdf",
        legend_position = (:left, :bottom),
        ykey = "entropy",
        ynorm = 1e-5,
        exponent_text = L"\times 10^{-8}", # TODO: check exponent
        xticks = [0, 3, 6, 9, 12],
        xlims = [0, 12],
        ylims = [-5, 0.25],
    )

    # Figure 5e
    plot_evolution(
        [
            "../results/20250520_steady_barotropic_instability_ec/N3M32",
            "../results/20250520_steady_barotropic_instability_es/N3M32",
        ],
        plot_name = "steady_barotropic_instability_entropy_evolution_N3M32.pdf",
        legend_position = (:left, :bottom),
        ykey = "entropy",
        ynorm = 1e-5,
        exponent_text = L"\times 10^{-8}", # TODO: check exponent
        xticks = [0, 3, 6, 9, 12],
        xlims = [0, 12],
        ylims = [-5, 0.25],
    )

    # Figure 5f
    plot_evolution(
        [
            "../results/20250520_steady_barotropic_instability_ec/N3M64",
            "../results/20250520_steady_barotropic_instability_es/N3M64",
        ],
        plot_name = "steady_barotropic_instability_entropy_evolution_N3M64.pdf",
        legend_position = (:left, :bottom),
        ykey = "entropy",
        ynorm = 1e-5,
        exponent_text = L"\times 10^{-8}",
        xticks = [0, 3, 6, 9, 12],
        xlims = [0, 12],
        ylims = [-5, 0.25],
    )

    # Figures 9a and 9b
    plot_timestep_study()
end

function plot_timestep_study()
    project_dir = joinpath(RESULTS_DIR, "20251219_steady_barotropic_instability_timestep_study")

    plot_timestep_study(project_dir)
end

function plot_timestep_study(
    project_dir::AbstractString;
    plots_dir = PLOTS_DIR,
    sort_rev = true,
    evolution_plot_name = "steady_barotropic_instability_timestep_entropy_evolution.pdf",
    convergence_plot_name = "steady_barotropic_instability_timestep_convergence.pdf",
)
    # Collect timestep run directories (expected suffix: _dt<value>)
    run_dirs = String[]
    dts = Float64[]
    for entry in readdir(project_dir; join = true)
        if isdir(entry) && occursin("_dt", entry)
            m = match(r"_dt(.+)$", entry)
            if m !== nothing
                dt = tryparse(Float64, m.captures[1])
                if dt !== nothing
                    push!(run_dirs, entry)
                    push!(dts, dt)
                end
            end
        end
    end

    if isempty(run_dirs)
        error(
            "plot_timestep_study: no timestep run directories found under $(project_dir). " *
            "Expected subdirectories matching '*_dt<value>'.",
        )
    end

    order = sortperm(dts; rev = sort_rev)
    run_dirs = run_dirs[order]
    dts = dts[order]
    labels = [LaTeXString("\$\\Delta t =\$" * @sprintf("%.6g\\,s", dt)) for dt in dts]
    colors = collect(Makie.cgrad(:blues, length(run_dirs); categorical = true))

    # Figure 9a: Entropy evolution over time (curves for different dt)
    plot_evolution(
        run_dirs;
        line_order = collect(1:length(run_dirs)),
        plots_dir = plots_dir,
        plot_name = evolution_plot_name,
        labels = labels,
        styles = fill(:solid, length(run_dirs)),
        colors = colors,
        ykey = "entropy",
        ylabel = LaTeXString("Normalized absolute entropy change"),
        relative = true,
        x_in_days = true,
        vlinepositions = Float64[],
        size = (500, 350),
        legend_position = (:left, :top),
        show_in_inset = collect(1:length(run_dirs)),
        xlims= [0, 12.0],
        xticks = [0, 3, 6, 9, 12],
        yscale = log10,
        plot_absolute = true,
        ylims = [1e-15, 1e-4],
    )

    # Figure 9b: Convergence plot at t = 12 days
    plot_convergence(
        [project_dir];
        plots_dir = plots_dir,
        plot_name = convergence_plot_name,
        file = "timestep_analysis.dat",
        select_cols = nothing,
        xkey = "dt",
        ykeys = ["entropy_error", "mass_error"],
        labels = [LaTeXString("Entropy"), LaTeXString("Mass")],
        styles = [:solid, :dash],
        colors = [1, 2],
        xlabel = L"Time step size $\Delta t$ (s)",
        ylabel = L"Normalized absolute change at $t = 12$ days",
        xscale = log10,
        yscale = log10,
        xticks = [10, 20, 40, 80, 160, 320],
        xlims = nothing,
        ylims = [1e-15, 1e-5],
        yticks = LogTicks(-15:-5),
        triangle_top = false,
        triangle_bottom = true,
        triangle_bottom_order = 5,
        legend_position = (:left, :top),
        size = (500, 350),
    )
end

function plot_rossby_haurwitz()
    # Figure 10a
    plot_evolution(
        [
            "../results/20250826_rossby_haurwitz_ec_N3/N3M16/",
            "../results/20250826_rossby_haurwitz_es_N3/N3M16/",
            "../results/20250826_rossby_haurwitz_standard_N3/N3M16/",
        ],
        line_order = [2, 1, 3],
        plot_name = "rossby_haurwitz_entropy_evolution_N3M16.pdf",
        labels = [LaTeXString("EC"), LaTeXString("ES"), LaTeXString("Standard")],
        styles = [:dash, :solid, :dot],
        ykey = "entropy",
        ynorm = 1e-6,
        exponent_text = L"\times 10^{-6}",
        xticks = [0, 7, 14, 21, 28],
        legend_position = (:right, :top),
        xlims = [0, 28],
        ylims = [-5, 5],
        vlinepositions = [18.71370, 26.03088], # crash times
        size = (500, 350),
    )

    # Figure 10b
    plot_evolution(
        [
            "../results/20250905_rossby_haurwitz_ec_N6/N6M8/",
            "../results/20250905_rossby_haurwitz_es_N6/N6M8/",
            "../results/20250905_rossby_haurwitz_standard_N6/N6M8/",
        ],
        line_order = [2, 1, 3],
        plot_name = "rossby_haurwitz_entropy_evolution_N6M8.pdf",
        labels = [LaTeXString("EC"), LaTeXString("ES"), LaTeXString("Standard")],
        styles = [:dash, :solid, :dot],
        ykey = "entropy",
        ynorm = 1e-6,
        exponent_text = L"\times 10^{-6}",
        xticks = [0, 7, 14, 21, 28],
        legend_position = (:right, :top),
        xlims = [0, 28],
        ylims = [-5, 5],
        vlinepositions = [13.39448, 18.01256], # crash times
        size = (500, 350),
    )

    # Figure 10c
    plot_evolution(
        [
            "../results/20250826_rossby_haurwitz_ec_N3/N3M16/",
            "../results/20250826_rossby_haurwitz_es_N3/N3M16/",
            "../results/20250826_rossby_haurwitz_standard_N3/N3M16/",
        ],
        line_order = [2, 1, 3],
        plot_name = "rossby_haurwitz_enstrophy_evolution_N3M16.pdf",
        labels = [LaTeXString("EC"), LaTeXString("ES"), LaTeXString("Standard")],
        styles = [:dash, :solid, :dot],
        ykey = "pot_enst",
        ylabel = LaTeXString("Normalized potential enstrophy change"),
        xticks = [0, 7, 14, 21, 28],
        legend_position = (:right, :top),
        xlims = [0, 28],
        ylims = [-2, 15],
        vlinepositions = [18.71370, 26.03088], # crash times
        size = (500, 350),
        inset_ylims = [-0.0011, 0.0001],
        inset_xlims = [0, 28],
        inset_xticks = [0, 7, 14, 21, 28],
        inset_yticks = [-0.001, 0],
        show_in_inset = [2],
    )

    # Figure 10d
    plot_evolution(
        [
            "../results/20250905_rossby_haurwitz_ec_N6/N6M8/",
            "../results/20250905_rossby_haurwitz_es_N6/N6M8/",
            "../results/20250905_rossby_haurwitz_standard_N6/N6M8/",
        ],
        line_order = [2, 1, 3],
        plot_name = "rossby_haurwitz_enstrophy_evolution_N6M8.pdf",
        labels = [LaTeXString("EC"), LaTeXString("ES"), LaTeXString("Standard")],
        styles = [:dash, :solid, :dot],
        ykey = "pot_enst",
        ylabel = LaTeXString("Normalized potential enstrophy change"),
        xticks = [0, 7, 14, 21, 28],
        legend_position = (:right, :top),
        xlims = [0, 28],
        ylims = [-2, 15],
        vlinepositions = [13.39448, 18.01256], # crash times
        size = (500, 350),
        inset_ylims = [-0.0011, 0.0001],
        inset_xlims = [0, 28],
        inset_xticks = [0, 7, 14, 21, 28],
        inset_yticks = [-0.001, 0],
        show_in_inset = [2],
    )
end