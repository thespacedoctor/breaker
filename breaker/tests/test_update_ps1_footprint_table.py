import os
import nose
import shutil
import yaml
from breaker import update_ps1_footprint_table, cl_utils
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


class test_update_ps1_footprint_table():

    def test_update_ps1_footprint_table_function(self):
        testObject = update_ps1_footprint_table(
            log,
            settings=settings,
            updateNed=False)
        testObject.get()

        # x-print-testpage-for-pessto-marshall-web-object

    # x-class-to-test-named-worker-function
