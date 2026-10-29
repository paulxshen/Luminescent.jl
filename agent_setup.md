# Lumi Agent Setup (Claude, ChatGPT)
### Claude
Experimental feature letting Claude use Lumi FDTD! Say for example "design 2 bus SOI220 ring resonator critically coupled Q>1000". Claude can optimize parameters and summarize results in a notebook. What can't Claude do?! Not yet production ready but something fun to try.

Any code the agent runs goes into notebooks in your Lumi Colab home at `~/projects/<name>/<notebook>.ipynb`. Related notebooks share a project, and each notebook's simulation files are in `~/projects/<name>/runs/<notebook>/`. When an agent (Claude, ChatGPT, DeepSeek, ...) updates a notebook, the previous version is kept in `~/projects/<name>/revisions/`.

On the free tier, jobs run one at a time.

claude.ai: Settings → Connectors → Add custom connector  
URL: https://api.luminescentai.com/mcp  
Sign in with your email at Lumi Colab, approve `lumi-mcp`, then return to Claude with the Lumi FDTD tools.
### Github Copilot
Like in VS Code, Lumi Colab has a Github Copilot sidebar that chats and edits notebook code!

### ChatGPT
Similar to Claude.

ChatGPT: Settings → Apps & Connectors → add a custom MCP server  
URL: https://api.luminescentai.com/mcp  
authentication: OAuth

### DeepSeek (untested)
DeepSeek Harness: add a remote (HTTP) MCP server  
URL: https://api.luminescentai.com/mcp  
header: `Authorization: Bearer <token>`, with a token from https://colab.luminescentai.com/hub/token  
DeepSeek's chat app doesn't support MCP.