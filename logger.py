import logging
import sys

# Logging
logger = logging.getLogger(__name__)
formatter = logging.Formatter("[%(filename)s:%(lineno)s - %(funcName)20s() - %(levelname)8s ] %(message)s")

console_handler: logging.StreamHandler = logging.StreamHandler(sys.stdout)
console_handler.setFormatter(formatter)
console_handler.setLevel(logging.INFO)

file_handler: logging.FileHandler = logging.FileHandler(
    filename="main.log",
    mode='w'
)
file_handler.setFormatter(formatter)
file_handler.setLevel(logging.DEBUG)

logger.addHandler(file_handler)
logger.addHandler(console_handler)
logger.setLevel(logging.INFO)

