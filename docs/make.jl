
using Documenter, Literate, ElementaryFluxModes

examples =
    sort(filter(x -> endswith(x, ".jl"), readdir(joinpath(@__DIR__, "src"), join = true)))

for example in examples
    Literate.markdown(
        example,
        joinpath(@__DIR__, "src"),
        repo_root_url = "https://github.com/HettieC/ElementaryFluxModes.jl/blob/main",
    )
end

example_mds = first.(splitext.(basename.(examples))) .* ".md"

withenv("COLUMNS" => 150) do
    makedocs(
        modules = [ElementaryFluxModes],
        clean = false,
        format = Documenter.HTML(
            ansicolor = true,
            canonical = "https://hettiec.github.io/ElementaryFluxModes.jl/stable/",
        ),
        sitename = "ElementaryFluxModes.jl",
        linkcheck = false,
        pages = ["README" => "index.md"; example_mds; "Reference" => "reference.md"],
    )
end

deploydocs(
    repo = "github.com/HettieC/ElementaryFluxModes.jl.git",
    target = "build",
    branch = "gh-pages",
    push_preview = false,
)
