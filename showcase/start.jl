# This script launches Pluto and opens a specific notebook file.
using Pluto
# using PlutoStyles

# Define the path to the notebook file.
# notebook_path = "/home/amartins/Documents/KUL/LineCableModels/showcase/showcase1.jl"

# println("Starting Pluto and opening notebook: ", notebook_path)
println("Starting Pluto, check your browser...")

try
	# Run Pluto, telling it which notebook to open
	Pluto.run(launch_browser = true)

	# Keep the script alive so the Pluto server doesn't shut down immediately.
	println("\nPluto server is running. Press Ctrl+C in this terminal to stop.")
	wait(Condition()) # Waits indefinitely until interrupted (Ctrl+C)

	println("\nAn error occurred while trying to run Pluto:")
	showerror(stdout, e)
	Base.show_backtrace(stdout, catch_backtrace())
finally
	println("\nPluto server stopped.")
end

println("Launcher script finished.")
