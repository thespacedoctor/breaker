import os
import nose
import unittest
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


class test_plot_wave_observational_timelines(unittest.TestCase):

    def test_plot_wave_observational_timelines_function(self):

        testObject = plot_wave_observational_timelines(
            log=log,
            settings=settings,
            plotType="timeline",
            gwid="G270580",
            projection="mercator",
            probabilityCut=False,
            telescope="atlas",
            timestamp=True)
        testObject.get()

    def test_generate_fits_image_map(self):
        plotter = plot_wave_observational_timelines(
            log=log,
            settings=settings,
            databaseConnRequired=False
        )

        # plotter.generate_fits_image_map(
        #     pathToProbMap=pathToInputDir + "/GW151226-bayestar.fits",
        #     gwid="G184098"
        # )

    def test_generate_pdf_image_map(self):
        plotter = plot_wave_observational_timelines(
            log=log,
            settings=settings,
            databaseConnRequired=False
        )

        plotter.generate_probability_plot(
            gwid="G184098",
            pathToProbMap=pathToInputDir + "/GW151226-bayestar.fits",
            fileFormats=["pdf", "png"],
            outputDirectory=pathToOutputDir,
            projection="mollweide",
            plotType="timeline",
            fitsImage=False,
            allSky=True,
            center=False
        )

        plotter.generate_probability_plot(
            gwid="G184098",
            pathToProbMap=pathToInputDir + "/GW151226-bayestar.fits",
            fileFormats=["pdf", "png"],
            outputDirectory=pathToOutputDir,
            projection="cartesian",
            plotType="timeline",
            fitsImage=False,
            allSky=True,
            center=False
        )

        plotter.generate_probability_plot(
            gwid="G184098",
            pathToProbMap=pathToInputDir + "/GW151226-bayestar.fits",
            fileFormats=["pdf", "png"],
            outputDirectory=pathToOutputDir,
            projection="mercator",
            plotType="timeline",
            fitsImage=False,
            allSky=True,
            center=False
        )

        # plotter.generate_probability_plot(
        #     gwid="G184098",
        #     pathToProbMap=pathToInputDir + "/GW151226-bayestar.fits",
        #     fileFormats=["pdf", "png"],
        #     outputDirectory=pathToOutputDir,
        #     projection="cartesian",
        #     plotType="timeline",
        #     fitsImage=False,
        #     allSky=True,
        #     center=False
        # )

        # x-print-testpage-for-pessto-marshall-web-object

    def test_generate_timeline_plot_maps(self):

        p = plot_wave_observational_timelines(
            log=log,
            settings=settings,
            gwid="G211117",
            plotType="timeline",
            allPlots=False,
            telescope="ps1",
            projection="mercator"
        )
        p.get()
