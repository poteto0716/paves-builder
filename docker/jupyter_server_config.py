"""Container defaults; keep Jupyter's random token unless explicitly configured."""
import os

c = get_config()  # noqa: F821 -- injected by Jupyter
c.ServerApp.ip = '0.0.0.0'
c.ServerApp.port = 8888
c.ServerApp.port_retries = 0
c.ServerApp.open_browser = False
c.ServerApp.root_dir = '/workspace'
c.ServerApp.default_url = '/lab/tree/workflows/adhesion/ovito_view.ipynb'
token = os.environ.get('PAVES_JUPYTER_TOKEN')
if token:
    c.IdentityProvider.token = token
