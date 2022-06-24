__version__ = "0.1.0"
import logging

log = logging.getLogger(__name__)
logging.basicConfig(
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    handlers=[logging.StreamHandler()],
)

from .ticfinder import TICFinder  # noqa
