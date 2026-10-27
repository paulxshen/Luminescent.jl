# Agent Setup
### Claude
Experimental feature letting Claude use Lumi FDTD! Say for example "design 2 bus SOI220 ring resonator critically coupled Q>1000". Claude can optimize parameters and summarize results in a notebook. What can't Claude do?! Not yet production ready but something fun to try.

Any code the gent runs go into notebooks in your Lumi Colab home at `~/projects/<name>/<name>.ipynb`. When an agent (Claude, ChatGPT, DeepSeek, ...) updates one, the previous version is kept in `~/projects/<name>/revisions/`.

claude.ai: Settings → Connectors → Add custom connector  
URL: https://api.luminescentai.com/mcp  
It should send you to GitHub (through the hub), then show the hub's "authorize lumi-mcp" page, then return to Claude with the 7 tools.
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
header: `Authorization: Bearer <token>`, with a token from https://hub.luminescentai.com/hub/token  
DeepSeek's chat app doesn't support MCP.