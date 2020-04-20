
using Pkg

println("Activating RetinalChaos files")
#Pkg.activate("RetinalChaos")
println("Developing RetinalChaos files")
#Pkg.develop("NeuroChaos")

function check_pkg(pkg::String)
    println("Checking package $pkg.")
    pkgs_installed = Pkg.installed()
    if haskey(pkgs_installed, pkg)
        println("Package $pkg is installed.")
    else
        println("Installing package $pkg.")
        Pkg.add(pkg)
        println("Package $pkg successfully installed.")
    end
    println("  ")
end

println("Checking RetinalChaos dependancies")
check_pkg("Dates")
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
