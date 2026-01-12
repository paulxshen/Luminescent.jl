# Installation

## Option 1: Luminescent Community cloud (free and public)
Run on our latest Nvidia Blackwell GPUs for free! Just [register](https://forms.gle/fP9wAkdJinT8t66w8) (allow 1 business day for processing) and install our python frontend client. Unlike before, no need to install the backend. See our simulation tutorials afterwards. Community cloud only supports simulations (not inverse design). Everything uploaded is publicly visible. If you want a private and secured environment, please upgrade to an enterprise cloud.


### Option 1A: Frontend client on Google Colab
`pip install luminescent`

Unlike before, we're only installing the python frontend client. This only takes a minute and the basic free CPU runtime is sufficient. See [Google Colab example](https://colab.research.google.com/drive/1g1sjXUTavK5ceZ-AeKyHQolpP2-RD0rd?usp=sharing)

### Option 1B: Frontend client on local machine
Due to a issue with `pymeshlab` dependency, we currently only support Python 3.12 . If that's your default python, do:  

`pip install luminescent`

If you installed multiple Python versions, you can also explicitly do:  

`path_to_python3.12 -m pip install luminescent` (replace `path_to_python3.12` with your actual python 3.12 path)

Some systems can take few minutes to initialize dependencies on first run. 

## Option 2: Luminescent Enterprise cloud (private and secure)
Similar to community cloud but private and secured for your organization. Comes with premium support and consulting. Pricing is a fixed annual fee plus $3 per GPU hour (similar to raw rates on GCP or AWS). Sign up via info@luminescentai.com

## Option 3: Local backend installation (advanced users)
Customers of Enterprise cloud are also entitled to local backend binaries upon request. 