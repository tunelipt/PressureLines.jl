using PressureLines
using Documenter

DocMeta.setdocmeta!(PressureLines, :DocTestSetup, :(using PressureLines); recursive=true)

makedocs(;
    modules=[PressureLines],
    authors="= <pjabardo@ipt.br> and contributors",
    sitename="PressureLines.jl",
    format=Documenter.HTML(;
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)
