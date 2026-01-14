from pyshockflow import Driver, Config
import sys
import traceback

try:
    config = Config('input.ini')
    driver = Driver(config)
    driver.solve()
    sys.exit(0)
except Exception as e:
    # Catch only real exceptions
    traceback.print_exc()
    sys.exit(1)