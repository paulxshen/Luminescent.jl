using Genie, Genie.Router, Genie.Assets

Genie.config.websockets_server = true # enable the websockets server

route("/") do
    Assets.channels_support() *
    """
    <script>
    window.parse_payload = function(payload) {
      console.log('Got this payload: ' + payload);
    }
    </script>
    """
end


channel("/____/echo") do
    @info "Received: $(params(:payload))"
end
up() # start the servers
# Genie.WebChannels sendMessageTo("____", "echo", "Hello!")
# Genie.WebChannels.broadcast("____", "Hey!")
Genie.WebChannels.connected_clients()
