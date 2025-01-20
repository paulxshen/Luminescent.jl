wget -nv "https://julialang-s3.julialang.org/bin/linux/x64/1.11/julia-1.11.2-linux-x86_64.tar.gz" -O /tmp/julia.tar.gz &> /dev/null
sudo tar -x -f /tmp/julia.tar.gz -C /usr/local --strip-components 1
julia -e 'using Pkg;pkg"add PackageCompiler"' &> /dev/null

julia -e 'using Pkg;pkg"up"'
pip install -U luminescent &> /dev/null

rm -rf Luminescent.jl
git clone https://paulxshen@github.com/paulxshen/Luminescent.jl &> /dev/null
cd Luminescent.jl

export PATH=~/LuminescentAI/bin:$PATH
julia "build/build.jl"