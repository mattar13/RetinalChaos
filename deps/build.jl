using Dates
println("[$(Dates.now())]: Building package RetinalChaos")
using Pkg

println("[$(Dates.now())]: Activating RetinalChaos files")
#Pkg.activate("RetinalChaos")
println("[$(Dates.now())]: Developing RetinalChaos files")
#Pkg.develop("NeuroChaos")

function check_pkg(pkg::String)
    println("[$(Dates.now())]: Checking package $pkg.")
    pkgs_installed = Pkg.installed()
    if haskey(pkgs_installed, pkg)
        println("[$(Dates.now())]: Package $pkg is installed.")
    else
        println("[$(Dates.now())]: Installing package $pkg.")
        Pkg.add(pkg)
        println("[$(Dates.now())]: Package $pkg successfully installed.")
    end
    println("  ")
end

println("[$(Dates.now())]: Checking RetinalChaos dependancies")
  
check_pkg("JSON2")
check_pkg("DifferentialEquations") 
check_pkg("ParameterizedFunctions")
check_pkg("LinearAlgebra") 
check_pkg("ForwardDiff") 
check_pkg("NLsolve")
check_pkg("Distributions")
check_pkg("Images") 
check_pkg("ImageSegmentation")
check_pkg("ProgressMeter")

#Pkg.resolve()
#Pkg.test()
#Pkg.activate()
