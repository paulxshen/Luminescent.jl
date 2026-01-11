# Installation

## Option 1: Luminescent Community Cluster (free and public)
Run on our latest Nvidia Blackwell GPUs for free! Just [register](https://forms.gle/fP9wAkdJinT8t66w8) (allow 1 business day for processing) and install our python frontend client. Unlike before, no need to install the backend. 

Due to a temporary issue with `pymeshlab` dependency, we currently only support Python 3.12 . If that's your default python, do:
`pip install luminescent`

If you installed multiple Python versions and don't know how to configure the default `python` using environmental variable PATH, you should explicitly do:
`path_to_python3.12 -m pip install luminescent` (replace `path_to_python3.12` with your actual python 3.12 path)

Some systems can take few minutes to initialize dependencies on first run. If your local machine doesn't work, you can use Google Colab which currently defaults to Python 3.12 . Unlike before, we're only installing the python frontend client. This only takes a minute and the basic free CPU runtime is sufficient. 

After installation, see our simulation tutorials. All simulations uploaded on the community cluster are publicly visible. If you want a private and secured environment, please upgrade to an enterprise cluster.

## Option 2: Luminescent Enterprise Cluster (private and secure)
Similar to community cluster but private and secured for your organization. Comes with premium support and consulting. Pricing is a fixed annual fee plus $3 per GPU hour (similar to raw rates on GCP or AWS). Sign up via info@luminescentai.com

## Option 3: Local Installation
Customers of Enterprise Cluster are also entitled to local installation binaries upon request. 