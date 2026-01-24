# Installation

## Community Tier: free public cloud
Run on our latest Nvidia Blackwell GPUs for free! Just [register](https://forms.gle/fP9wAkdJinT8t66w8) (allow 1 business day for processing) and install our python frontend client. Unlike before, no need to install the backend. See our simulation tutorials afterwards. Community cloud only supports simulations (not inverse design). Everything uploaded is publicly visible. 

### Option A: Frontend client on Google Colab
`pip install luminescent`

Unlike before, we're only installing the python frontend client. This only takes a minute and the basic free CPU runtime is sufficient. See [Google Colab example](https://colab.research.google.com/drive/1g1sjXUTavK5ceZ-AeKyHQolpP2-RD0rd?usp=sharing)

### Option B: Frontend client on local machine
Due to a issue with `pymeshlab` dependency, we currently only support Python 3.12 . If that's your default python, do:  

`pip install luminescent`

If you installed multiple Python versions, you can also explicitly do:  

`path_to_python3.12 -m pip install luminescent`

Some systems can take few minutes to initialize dependencies on first run. 

## Enterprise Tier: private cloud or local install
Pay once for Enterprise Tier and run forever! You'll get backend binaries for running simulations on your own machine. Alternatively, if you love the hassle free cloud experience, you can now run jobs privately and securely on our Enterprise Cloud with a small GPU hour surcharge equivalent to raw cost from AWS or GCP. Sign up via info@luminescentai.com
