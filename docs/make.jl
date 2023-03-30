using Documenter
using OSQPSolver

makedocs(
    sitename = "OSQPSolver",
    format = Documenter.HTML(),
    modules = [OSQPSolver]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
