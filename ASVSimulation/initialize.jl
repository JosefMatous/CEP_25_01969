# Initialize the ASV Simulation project by resolving local dependencies and precompiling the packages.

import Pkg
import TOML

# Find all local packages in the working directory
isdir_local(fname::String) = isdir(joinpath(@__DIR__, fname))
packages = filter(isdir_local, readdir(@__DIR__, join=false))

# Find local dependencies in the `Project.toml` file
local_deps = Dict{String, Set{String}}()
local_filter(name::String) = name in packages
for package in packages
    project_path = joinpath(@__DIR__, package, "Project.toml")
    if isfile(project_path)
        project = TOML.parsefile(project_path)
        deps = get(project, "deps", Dict{String, Any}())
        local_deps[package] = filter(local_filter, keys(deps))
    else
        error("Project.toml not found for package $package at path $project_path")
    end
end

# Iteratively resolve local dependencies
resolved = Set{String}()
is_resolved(package) = package in resolved
is_resolvable(package) = all(is_resolved, local_deps[package])

for _ in 1:length(packages)
    if all(is_resolved, packages)
        break
    end
    for package in packages
        if !is_resolved(package) && is_resolvable(package)
            Pkg.activate(joinpath(@__DIR__, package))
            for dep in local_deps[package]
                Pkg.develop(path=joinpath(@__DIR__, dep))
            end
            Pkg.instantiate()
            Pkg.precompile()
            push!(resolved, package)
        end
    end
end
if !all(is_resolved, packages)
    error("Could not resolve all local dependencies. Unresolved packages: $(setdiff(packages, resolved))")
end
     