using Documenter
using MinimalistOSQP

makedocs(
    sitename = "MinimalistOSQP",
    format = Documenter.HTML(),
    modules = [MinimalistOSQP]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
