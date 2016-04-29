import os
import nose
import shutil
import yaml
from breaker import plot_multi_panel_alternate_map_comparison, cl_utils
from breaker.utKit import utKit

from fundamentals import tools, times

su = tools(
    arguments={},
    docString=__doc__,
    logLevel="DEBUG",
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


class test_plot_multi_panel_alternate_map_comparison():

    def test_plot_multi_panel_alternate_map_comparison_function(self):
        kwargs = {}
        kwargs["log"] = log
        kwargs["settings"] = settings
        # xt-kwarg_key_and_value

        testObject = plot_multi_panel_alternate_map_comparison(**kwargs)
        testObject.get()

        # x-print-testpage-for-pessto-marshall-web-object

    # x-class-to-test-named-worker-function
