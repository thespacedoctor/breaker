#!/usr/local/bin/python
# encoding: utf-8
"""
*Generate a catalogue of fake transient sources for GW*
"""
################# GLOBAL IMPORTS ####################
import sys
import os
os.environ['TERM'] = 'vt100'
from HMpTy import htm
from fundamentals import tools, times
from fundamentals.mysql import readquery, table_exists, writequery
from HMpTy.mysql import add_htm_ids_to_mysql_database_table
import numpy as np
import math
import csv


class generate_faker_catalogue():
    """
    *Generate a text file with fake transients locations and magnitudes given a PS1 exposure ID and a gravitational wave ID*

    **Key Arguments:**
        - ``log`` -- logger
        - ``settings`` -- the settings dictionary
        - ``gwid`` -- the gravitational wave ID
        - ``ps1ExpId`` -- the PS1 footprint ID

    **Usage:**

        To populate the associated galaxies with fake sources and writing faker results to file run the code:

        .. code-block:: python 

            from breaker import fakers
            fakecat = fakers.generate_faker_catalogue(
                log=log,
                settings=settings,
                ps1ExpId=73283973,
                gwid="G211117"
            ).get() 

        The location of the output faker file is determined from the ``output directory`` setting in the settings file:

        .. code-block:: yaml 

            output directory: "/Users/Dave/Desktop"
    """

    def __init__(
            self,
            log,
            gwid,
            settings=False,
            ps1ExpId=False
    ):
        self.log = log
        log.debug("instansiating a new 'generate_faker_catalogue' object")
        self.settings = settings
        self.gwid = gwid
        self.ps1ExpId = ps1ExpId
        # xt-self-arg-tmpx

        # Initial Actions
        # CONNECT TO THE VARIOUS DATABASES REQUIRED
        from breaker import database
        db = database(
            log=self.log,
            settings=self.settings
        )
        self.ligo_virgo_wavesDbConn, self.ps1gwDbConn, self.cataloguesDbConn = db.get()

        return None

    def get(self):
        """
        *Generate the fake transient source file*

        **Return:**
            - None
        """
        self.log.info('starting the ``get`` method')

        self.generate_faker_locations_within_galaxies_found_in_survey_footprint()
        self.write_fakers_to_file()

        self.log.info('completed the ``get`` method')
        return None

    def write_fakers_to_file(
            self):
        """
        *Write out the associated fake transient sources to file, alongside 15% extra transients added blindly throughout the survey field*

        .. warning::

            The code will fail if sherlock has not yet been run against the survey areas for the give GW ID.

        **Usage:**

            .. code-block:: python 

                from breaker import fakers
                fakecat = fakers.generate_faker_catalogue(
                    log=log,
                    settings=settings,
                    ps1ExpId=73283973,
                    gwid="G211117"
                ).write_fakers_to_file() 
        """
        self.log.info('starting the ``write_fakers_to_file`` method')

        mesh16 = htm.HTM(16)
        theseArrays = []
        radius = 1.45
        ra = []
        dec = []

        gwid = self.gwid
        ps1ExpId = self.ps1ExpId

        if gwid:
            gwidClause = 'and gw_id = "%(gwid)s"' % locals()
        else:
            gwidClause = ""

        # SELECT THE BORE-SIGHT FOR THE PS1 EXPOSURE
        sqlQuery = u"""
              select raDeg, decDeg, gw_id from ps1_pointings where 1=1 %(gwidClause)s  and ps1_exp_id = %(ps1ExpId)s
        """ % locals()
        rows = readquery(
            log=self.log,
            sqlQuery=sqlQuery,
            dbConn=self.ligo_virgo_wavesDbConn
        )
        raFP = rows[0]["raDeg"]
        decFP = rows[0]["decDeg"]
        gwid = rows[0]["gw_id"]

        # GENERATE AN ARRAY OF HTM IDs FOR THE FOV
        thisArray = mesh16.intersect(
            raFP, decFP, radius, inclusive=True)
        hmax = thisArray.max()
        hmin = thisArray.min()
        ratio = float(hmax - hmin + 1) / float(thisArray.size)
        if ratio < 100 or thisArray.size > 2000:
            htmWhereClause = "where htm16ID between %(hmin)s and %(hmax)s" % locals(
            )
        else:
            s = StringIO()
            np.savetxt(s, thisArray, fmt='%d', newline=",")
            thesHtmIds = s.getvalue()[:-1]
            htmWhereClause = "where htm16ID in (%(thesHtmIds)s)" % locals()

        # TEST TABLE EXIST
        thisTableExists = table_exists(
            dbConn=self.ligo_virgo_wavesDbConn,
            log=self.log,
            dbTableName="tcs_%(gwid)s_catalogued_sources" % locals()
        )
        if not thisTableExists:
            self.log.error(
                'tcs_%(gwid)s_catalogued_sources table does not exist in LV database - need to run sherlock to create this table' % locals())
            return

        # ADD THE HTMIDs TO THE RELEVANT CATALOGUE SOURCES TABLE BEFORE
        # QUERYING
        add_htm_ids_to_mysql_database_table(
            raColName="catalogue_object_ra",
            declColName="catalogue_object_dec",
            tableName="tcs_%(gwid)s_catalogued_sources" % locals(),
            dbConn=self.ligo_virgo_wavesDbConn,
            log=self.log,
            primaryIdColumnName="id",
            batchSize=50000
        )

        # FINALLY SELECT THE FAKER DETAILS FROM FOV
        sqlQuery = """select faker_ra, faker_dec, z, catalogue_object_id, catalogue_object_subtype, 2mass_k_mag, 2mass_k_mag_error from tcs_%(gwid)s_catalogued_sources %(htmWhereClause)s and z is not null and z < 0.15 and (z_quality is null or z_quality not like 'PHOT%%') and (catalogue_object_subtype is null or catalogue_object_subtype not like "%%*%%") and major_axis_arcsec is not null """ % locals(
        )
        rows = readquery(
            log=self.log,
            sqlQuery=sqlQuery,
            dbConn=self.ligo_virgo_wavesDbConn
        )

        savedRows = []
        raList = []
        decList = []
        for row in rows:
            raList.append(row["faker_ra"])
            decList.append(row["faker_dec"])

        tRa = np.array([raFP])
        tDec = np.array([decFP])
        raList = np.array(raList)
        decList = np.array(decList)
        indexList1, indexList2, separation = mesh16.match(
            tRa, tDec, raList, decList, radius, maxmatch=0)
        for i in range(indexList1.size):
            savedRows.append(rows[indexList2[i]])

        # ADD IN AN EXTRA 15% FOR FAKERS FOUND IN THE FIELD BUT NOT NEAR GALAXY
        extraRows = []
        extra15 = int(len(savedRows) * 0.176)
        pi = (4 * math.atan(1.0))
        DEG_TO_RAD_FACTOR = pi / 180.0
        RAD_TO_DEG_FACTOR = 180.0 / pi
        for i in range(extra15):
            theta = 360. * np.random.rand()
            sep = (np.random.rand() * radius**2.)**(0.5)
            dec2 = decFP - sep * math.cos(theta * DEG_TO_RAD_FACTOR)
            ra2 = raFP - sep * math.sin(theta * DEG_TO_RAD_FACTOR) / math.cos(
                (2 * decFP - sep * math.cos(theta * DEG_TO_RAD_FACTOR)) * DEG_TO_RAD_FACTOR / 2)
            extraRows.append({
                "faker_ra": ra2,
                "faker_dec": dec2,
                "z": 0.00000000,
                "catalogue_object_id": "ORPHAN",
                "2mass_k_mag": 0.00000000,
                "2mass_k_mag_error": 0.00000000
            })

        rows = savedRows + extraRows
        e = 2.7182818284590452353602875
        writeRows = []
        # CALCULATE THE FAKER RANDOM MAGNITUDES
        for row in rows:
            mag = (np.random.rand() * 3.) + 17.
            f0 = 10**(-(mag + 48.6) / (2.5))
            ft = 0.25 * f0 * (e**2) * ((3.22 / 1.18)**2) * (e**(-3.22 / 1.18))
            magt = (-2.5) * (math.log10(ft)) - 48.6
            if not row["2mass_k_mag"]:
                row["2mass_k_mag"] = 0.00
            if not row["2mass_k_mag_error"]:
                row["2mass_k_mag_error"] = 0.00
            writeRows.append({
                "ra": row["faker_ra"],
                "dec": row["faker_dec"],
                "imag": magt,
                "z": row["z"],
                "catalogue_object_id": row["catalogue_object_id"],
                "2mass_k_mag": row["2mass_k_mag"],
                "2mass_k_mag_error": row["2mass_k_mag_error"]
            })

        # WRITE THE RESULT TO FILE
        import codecs
        from datetime import datetime, date, time
        now = datetime.now()
        now = now.strftime("%Y%m%dt%H%M%S")

        # RECURSIVELY CREATE MISSING DIRECTORIES - PATH FROM SETTINGS FILE
        if not os.path.exists(self.settings["output directory"] + "/fakers"):
            os.makedirs(self.settings["output directory"] + "/fakers")

        writeFile = codecs.open(
            self.settings["output directory"] + "/fakers/G184098_fake_sources_%(now)s_trimmed.csv" % locals(), encoding='utf-8', mode='w')
        writeFile2 = codecs.open(
            self.settings["output directory"] + "/fakers/G184098_fake_sources_%(now)s_complete.csv" % locals(), encoding='utf-8', mode='w')
        writeFile.write("index,ra,dec,i-mag\n")
        writeFile2.write(
            "index,ra,dec,i-mag,redshift,galaxy-id,2mass-k-mag,2mass-k-mag-error\n")
        for r, d in enumerate(writeRows):
            d["r"] = r + 1
            one = "%(r)04d,%(ra)6.5f,%(dec)6.5f,%(imag)4.2f\n" % d
            writeFile.write(one.replace(
                ",0.00,", ",null,").replace(",0.00\n", ",null\n"))
            two = '%(r)04d,%(ra)6.5f,%(dec)6.5f,%(imag)4.2f,%(z)4.3f,"%(catalogue_object_id)s",%(2mass_k_mag)4.2f,%(2mass_k_mag_error)4.2f\n' % d
            writeFile2.write(two.replace(
                ",0.00,", ",null,").replace(",0.00\n", ",null\n"))
        writeFile.close()
        writeFile2.close()

        self.log.info('completed the ``write_fakers_to_file`` method')
        return None

    def generate_faker_locations_within_galaxies_found_in_survey_footprint(
            self):
        """
        *For galaxies matched within the PS1 survey footprints, add faker sources at a random location within the galaxy*

        .. note::

            the galaxies need to have a semi-major axis measurement before a faker is added

        **Return:**
            - None

        **Usage:**

            .. code-block:: python 

                # ADD FAKERS TO GALAXIES
               from breaker import fakers
                fakecat = fakers.generate_faker_catalogue(
                    log=log,
                    settings=settings,
                    ps1ExpId=73283973,
                    gwid="G211117"
                ).generate_faker_locations_within_galaxies_found_in_survey_footprint() 
        """
        self.log.info(
            'starting the ``generate_faker_locations_within_galaxies_found_in_survey_footprint`` method')

        import math
        pi = (4 * math.atan(1.0))
        DEG_TO_RAD_FACTOR = pi / 180.0
        RAD_TO_DEG_FACTOR = 180.0 / pi

        # ITERATE THROUGH THE CATALOGUED SOURCES DATABASE TABLES
        for gwid in self.settings["gravitational waves"]:
            tableName = "tcs_%(gwid)s_catalogued_sources" % locals()
            thisTableExists = table_exists(
                dbConn=self.ligo_virgo_wavesDbConn,
                log=self.log,
                dbTableName=tableName
            )
            if thisTableExists:
                # SELECT GALAXY ROWS WITH SEMI-MAJOR AXIS MEASUREMENT BUT NO
                # FAKER COORDINATES
                sqlQuery = u"""
                    select id, catalogue_object_ra, catalogue_object_dec, major_axis_arcsec from %(tableName)s where  major_axis_arcsec is not null and faker_ra is null
                """ % locals()
                rows = readquery(
                    log=self.log,
                    sqlQuery=sqlQuery,
                    dbConn=self.ligo_virgo_wavesDbConn
                )
                count = 0
                values = []

                if len(rows) == 0:
                    print "Faker sources have been added to all galaxies found within the PS1 footprint"

                # GENERATE A RANDOM FAKER LOCATION AND ADD TO THE DATABASE -
                # ADD INT BATCHES OF 5000
                for row in rows:
                    count += 1
                    thisId = int(row["id"])
                    major_axis_arcsec = row["major_axis_arcsec"]
                    sep = (major_axis_arcsec /
                           3660.) * np.random.rand()
                    arcsecSep = sep * 3600.
                    print np.random.rand()
                    theta = 360. * np.random.rand()

                    ra1 = row["catalogue_object_ra"]
                    dec1 = row["catalogue_object_dec"]

                    dec2 = dec1 - sep * math.cos(theta * DEG_TO_RAD_FACTOR)
                    ra2 = ra1 - sep * math.sin(theta * DEG_TO_RAD_FACTOR) / math.cos(
                        (2 * dec1 - sep * math.cos(theta * DEG_TO_RAD_FACTOR)) * DEG_TO_RAD_FACTOR / 2)

                    values.append(str((thisId, ra2, dec2)))
                    if count == 5000:
                        values = ",".join(values)
                        sqlQuery = u"""
                            INSERT INTO %(tableName)s (id,faker_ra,faker_dec) VALUES %(values)s
                            ON DUPLICATE KEY UPDATE faker_ra=VALUES(faker_ra),faker_dec=VALUES(faker_dec);
                        """ % locals()
                        writequery(
                            log=self.log,
                            sqlQuery=sqlQuery,
                            dbConn=self.ligo_virgo_wavesDbConn
                        )
                        count = 0
                        values = []

                        print "Another 5000 faker sources have been added to galaxies found within the PS1 footprint"

                values = ",".join(values)

                if len(values):
                    sqlQuery = u"""
                        INSERT INTO %(tableName)s (id,faker_ra,faker_dec) VALUES %(values)s
                        ON DUPLICATE KEY UPDATE faker_ra=VALUES(faker_ra),faker_dec=VALUES(faker_dec);
                    """ % locals()
                    writequery(
                        log=self.log,
                        sqlQuery=sqlQuery,
                        dbConn=self.ligo_virgo_wavesDbConn
                    )

                    count = len(values)
                    print "Another %(count)s faker sources have been added to galaxies found within the PS1 footprint" % locals()

        self.log.info(
            'completed the ``generate_faker_locations_within_galaxies_found_in_survey_footprint`` method')
        return None
