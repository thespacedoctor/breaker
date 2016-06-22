import os
import nose
import shutil
import yaml
from breaker import cl_utils
from breaker.utKit import utKit
from breaker.plots import plot_wave_observational_timelines

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


class test_plot_wave_observational_timelines():

    def test_plot_wave_observational_timelines_function(self):

        testObject = plot_wave_observational_timelines(
            log=log,
            settings=settings,
            plotType="timeline",
            gwid="G184098",
            projection="tan",
            probabilityCut=False)
        testObject.get()

        # x-print-testpage-for-pessto-marshall-web-object

    # x-class-to-test-named-worker-function
