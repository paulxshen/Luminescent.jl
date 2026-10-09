for x=["bend", "ring", "TFLN_coupler", "GC", "GC_topopt"]
    cp("../moon/luminescent/$x.ipynb", "$x.ipynb"; force=true)
end