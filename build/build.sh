wget -nv "https://julialang-s3.julialang.org/bin/linux/x64/1.11/julia-1.11.2-linux-x86_64.tar.gz" -O /tmp/julia.tar.gz &> /dev/null
sudo tar -x -f /tmp/julia.tar.gz -C /usr/local --strip-components 1
julia -e 'using Pkg;pkg"add PackageCompiler;up"' &> /dev/null

julia -e 'using Pkg;pkg"up"'
pip install -U luminescent &> /dev/null

rm -rf Luminescent.jl
git clone https://paulxshen@github.com/paulxshen/Luminescent.jl 

cd Luminescent.jl
julia "build/build.jl"
cd ..
alias julia='julia -J/usr/local/
# ../Luminescent/bin/Luminescent
# export PATH=~/LuminescentAI/bin:$PATH
export PATH=Luminescent/bin:$PATH
lumijulia -e 'using Pkg;pkg"add CUDA"'
tar czf Luminescent.tar.gz Luminescent
# split -b 1500M Luminescent.tar.gz "Luminescent.tar.gz.part"