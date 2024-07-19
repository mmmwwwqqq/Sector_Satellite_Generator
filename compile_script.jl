using PackageCompiler

# Path to your project directory
project_dir = "D:/Desktop/АЦЦКИЙ ПТСР V2.0/_code"

# Path to your compile script if needed
compile_script_path = joinpath(project_dir, "compile_script.jl")

# Create the application
create_app(project_dir, "MyApp"; script=compile_script_path)
