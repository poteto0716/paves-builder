"""Check that Jupyter responds, without putting its authentication token in argv."""
from urllib.request import urlopen

with urlopen('http://127.0.0.1:8888/login', timeout=5) as response:
    if response.status != 200:
        raise SystemExit(1)
