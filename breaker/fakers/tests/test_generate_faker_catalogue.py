import os
import nose
import shutil
import yaml
from breaker import cl_utils
from breaker.fakers import generate_faker_catalogue
from breaker.utKit import utKit

from fundamentals import tools, times

su = tools(
    arguments={},
    docString=__doc__,
    logLevel="WARNING",
    options_first=False,
    projectName="breaker"
)
arguments, settings, log, dbConn = su.setup()

# load settings
stream = file(
    "/Users/Dave/.config/breaker/breaker.yaml", 'r')
settings = yaml.load(stream)
stream.close()

# SETUP AND TEARDOWN FIXTURE FUNCTIONS FOR THE ENTIRE MODULE
moduleDirectory = os.path.dirname(__file__)
utKit = utKit(moduleDirectory)
log, dbConn, pathToInputDir, pathToOutputDir = utKit.setupModule()
utKit.tearDownModule()


class test_generate_faker_catalogue():

    def test_generate_faker_catalogue_function(self):
        kwargs = {}
        kwargs["log"] = log
        kwargs["settings"] = settings
        kwargs["gwid"] = "G184098"
        kwargs["ps1ExpId"] = 65925995
        # xt-kwarg_key_and_value

        testObject = generate_faker_catalogue(**kwargs)
        testObject.get()

        # x-print-testpage-for-pessto-marshall-web-object

    # x-class-to-test-named-worker-function
